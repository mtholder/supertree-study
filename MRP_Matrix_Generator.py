#!/usr/bin/env python
import sys
import os
import subprocess
from dendropy.treesplit import encode_splits,split_to_list
from cStringIO import StringIO
from dendropy.utility.messaging import get_logger
from dendropy.treecalc import fitch_down_pass, fitch_up_pass
from dendropy import DataSet
from dendropy.utility.error import DataParseError
from dendropy.utility.textutils import escape_nexus_token
_DEBUGGING = True
_LOG = get_logger('geodispersal')
verbose = False
AREA_NAME_LIST = []
col_width = 17

def warn(msg):
    _LOG.warn(msg)
LAST_COMMAND = ''
def write_as_nexus(stream, patterns, label):
    global LAST_COMMAND
    stream.write("\n[!%s ]\n" % label)
    p = patterns[0]
    num_chars = len(patterns)
    num_areas = len(p)
    if num_areas < len(AREA_NAME_LIST):
        warn('%d labels were found in the labels file, but only %d areas were read in the input NEXUS files' % (
                    len(AREA_NAME_LIST),
                    num_areas))
    elif num_areas > len(AREA_NAME_LIST):
        warn('Only %d labels were found in the labels file, but %d areas were read in the input NEXUS files' % (
                    len(AREA_NAME_LIST),
                    num_areas))
    stream.write("""Begin Data;
    Dimensions ntax = %d nchar = %d;
    Format datatype=standard symbols="012" ;
    Matrix \n""" % (num_areas, num_chars))
    for area_ind in range(num_areas):
        if AREA_NAME_LIST:
            try:
                name = AREA_NAME_LIST[area_ind]
            except:
                name = 'unlabelled area ' + str(1 + area_ind)
        else:
            name = "area%d" % area_ind
        name = escape_nexus_token(name)
        padding = ' ' * (col_width - len(name))
        stream.write("%s%s" % (name, padding))
        for p in patterns:
            stream.write(" %s" % str(p[area_ind]))
        stream.write("\n")
    if num_areas > 11:
        search = "HSearch"
    else:
        search = "BAndB"
    stream.write(""";
End;
Begin Paup;
    typeset * o = ord : 1-. ;
    %s ;
    SaveTrees from =1 to=100 file = %s.tre;
    %s
End;
""" % (search, label, LAST_COMMAND))
    stream.flush()


def vicariance_patterns(node_list, num_areas):
    '''Returns a transposed (characters as rows) matrix of 0, 1, 2 values for the
    vicariance matrix'''

    patterns = []
    root = node_list[0]
    assert(root.parent_node is None)
    curr_pat = [0]*num_areas
    for i in root.state_sets[0]:
        curr_pat[i] = 1
    patterns.append(curr_pat)
    areas_seen = set()
    count = 0
    for node in node_list[1:]:
        curr_pat = [0]*num_areas
        p = node.parent_node
        node.biogeo_number = count
        assert(p)
        par_area = p.state_sets[0]
        child_area = node.state_sets[0]
        areas_seen.update(par_area)
        areas_seen.update(child_area)

        diff = par_area - child_area
        _LOG.debug("count = %d Par = %s       Des = %s    diff = %s" %(count, str(par_area), str(child_area), str(diff)))
        if diff:
            # range contraction...
            # parent existed in areas that the descendant did not
            twos = set.intersection(child_area, par_area)
            ones = set.symmetric_difference(child_area, par_area)
            for i in twos:
                curr_pat[i] = 2
        else:
            ones = set.union(child_area, par_area)
        for i in ones:
            curr_pat[i] = 1
        _LOG.debug("curr_pat = %s" %(str(curr_pat)))
        patterns.append(curr_pat)
        count += 1

    if ABSENT_FOR_ALLTAXA_CODE != 0:
        replace_unseen_areas(patterns, num_areas, areas_seen)
    return patterns

def replace_unseen_areas(patterns, num_areas, areas_seen):
    for area_ind in range(num_areas):
        if area_ind not in areas_seen:
            for pattern in patterns:
                pattern[area_ind] = ABSENT_FOR_ALLTAXA_CODE

def dispersal_patterns(node_list, num_areas):
    '''Returns a transpose (characters as rows) matrix of 0, 1, 2 values for the
    dispersal matrix'''
    patterns = []
    root = node_list[0]
    assert(root.parent_node is None)
    curr_pat = [0]*num_areas
    for i in root.state_sets[0]:
        curr_pat[i] = 1
    patterns.append(curr_pat)

    areas_seen = set()

    for node in node_list[1:]:
        p = node.parent_node
        assert(p)
        par_area = p.state_sets[0]
        child_area = node.state_sets[0]

        areas_seen.update(par_area)
        areas_seen.update(child_area)

        twos = child_area - par_area
        curr_pat = [0]*num_areas
        for i in twos:
            curr_pat[i] = 2
        for i in par_area: #BRUCE's EMAIL mod...  par_area instead of: set.intersection(child_area, par_area):        for i in set.intersection(child_area, par_area):
            curr_pat[i] = 1
        _LOG.debug("dispersal: count = %d Par = %s, Des = %s, twos = %s, pattern= %s" % 
                    (node.biogeo_number, str(par_area), str(child_area), str(twos), str(curr_pat)))
        patterns.append(curr_pat)

    if ABSENT_FOR_ALLTAXA_CODE != 0:
        replace_unseen_areas(patterns, num_areas, areas_seen)
    return patterns

ABSENT_FOR_ALLTAXA_CODE = 0

def add_to_full(mat, full_mat, mat_el_len_list):
    prev_n_char = len(full_mat)
    prev_n_areas = prev_n_char > 0 and len(full_mat[0]) or 0
    n_additional_char = len(mat)
    assert(n_additional_char > 0)
    n_areas_mat = len(mat[0])
    if n_areas_mat > prev_n_areas:  
        for el in full_mat:
            el.extend([ABSENT_FOR_ALLTAXA_CODE]*(n_areas_mat - prev_n_areas))
    elif n_areas_mat < prev_n_areas:  
        for el in mat:
            el.extend([ABSENT_FOR_ALLTAXA_CODE]*(prev_n_areas - n_areas_mat))

    mat_el_len_list.append(n_additional_char)
    full_mat.extend(mat)

def check_file_overwrite(fn):
    if os.path.exists(fn):
        c = raw_input(fn  + ' exists. Are you sure that you would like to over-write it? (y/n)\n')
        if c.lower() != 'y':
            sys.exit('Did not get "y" in response to prompt to overwrite file. Exiting...')

if __name__ == '__main__':
    output = sys.stdout
    fo = open(sys.argv[1], "rU")
    dataset = DataSet()
    try:
        dataset.read(stream=fo, schema="Newick")
    except DataParseError as dfe:
        raise ValueError(str(dfe))
    if len(dataset.taxon_sets) != 1:
        raise ValueError("Expecting one set of taxa in %s" % f)
    if len(dataset.tree_lists) != 1:
        raise ValueError("Expecting one tree in %s" % f)
    taxon_set = dataset.taxon_sets[0]
    tree_list = dataset.tree_lists[0]
    #print len(tree_list)
    number_of_taxon = len(taxon_set)
    branch_counter = 0
    code_list = [StringIO() for i in taxon_set]
    for tree in tree_list:
        encode_splits(tree)
        tree_mask = tree.seed_node.edge.split_bitmask
        assert tree_mask is not None
        tree_tax = set(split_to_list(tree_mask))
        print tree_tax
        split_list = []
        for node in tree.postorder_internal_node_iter():
            if node.parent_node is not None:
                branch_counter +=1
                split_set = set(split_to_list(node.edge.split_bitmask))
                split_list.append(split_set)
        for i,stream in enumerate(code_list):
            if i in tree_tax:
                for split in split_list:
                    if i in split:
                        stream.write('1')
                    else:
                        stream.write('0')
            else:
                stream.write('?'*len(split_list))
                        
            
    #print branch_counter
    output.write("""#NEXUS
    Begin Data;
    Dimensions ntax = %d nchar = %d;
    Format datatype=standard symbols="01" Missing = ?;
    Matrix \n""" % (number_of_taxon, branch_counter))
    
    escaped_names = []
    max_name_length = 0
    for t in taxon_set:
        escaped_names.append(escape_nexus_token(t.label))
        max_name_length = max(max_name_length, len(escaped_names[-1]))
    fmt = '%%-%ds %%s\n' %(max_name_length)
    for i,stream in enumerate(code_list):
        output.write(fmt %(escaped_names[i],stream.getvalue()))

    output.write(""";
    END;
    """)
        

