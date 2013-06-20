import sys
import dendropy
from dendropy.treesplit import encode_splits,split_to_list
from dendropy.utility.error import DataParseError
import StringIO
class Taxon(object):
    def __init__(self,name,parent):
        self.name = name
        self.parent = parent
        self.children = []
        
    def write_as_newick(self, out):
        #if self.name == 'Primates':
        #    out = sys.stderr
        if self.children:
            out.write('(')
            for n,c in enumerate(self.children):
                if n > 0:
                    out.write(',')
                c.write_as_newick(out)
            out.write(')')
        out.write("'" + self.name + "'")

def parse_ott(filename = 'ott2-taxo-first500.txt',splitchar ='\t|\t'):
    '''
    function to be used..
    filename is based on user input, but is also
    given by a default value.
    
    splitchar is given a default value; however,
    this can be changed via function input.
    '''
    root = None
    keyToObject = {}
    with open(filename, 'rU') as inp_file:
        inp = iter(inp_file)
        next(inp)
        for row in inp:
            row_list = row.split(splitchar)
            t_id = int(row_list[0])
            name = row_list[2]
            ottoid = t_id
            mod_name = name + '@' + str(ottoid)
            mod_name = mod_name.replace("'","''")
            #if name is not "Microcebus sp. d'Ambre":
            #name.replace("'","''")
            try:
                parent_id = int(row_list[1])
                parent = keyToObject.get(parent_id)
            except:
                parent = None      
            taxon = Taxon(mod_name, parent)
            keyToObject[t_id] = taxon
            if root is None:
                root = taxon
            if parent is not None:
                #print('parent of ' + t_id + ' is ' + parent.name)
                parent.children.append(taxon)
            else:
                pass
                #print('No parent for ' + name + 'looked up ' + str(parent_id))
                #print(str(keyToObject))          
    return root



def modded_parse_ott(filename,splitchar = '\t|\t'):
    '''
This function is a modified version of parse_ott.
The filename is based on user input as is splitchar.
However splitchar is given a default value, which can be
overridden by user input.

This function is just to return the species primates and
their evolutionary offspring (for now at least..)
'''
    root = None
    keyToObject = {}
    primate_id_set = set()
    with open(filename, 'rU') as inp_file:
        inp = iter(inp_file)
        next(inp)
        for row in inp:
            row_list = row.split(splitchar)
            t_id = int(row_list[0])
            name = row_list[2]
            parent_id = None
            try:
                parent_id = int(row_list[1])
            except:
                pass
            if parent_id is not None:
                if name == 'Primates' or parent_id in primate_id_set:
                    sys.stdout.write(row)
                    primate_id_set.add(t_id)


def life_ott(fo):
    output = sys.stdout
    dataset = dendropy.DataSet()
    try:
        dataset.read(stream=fo, schema="Newick")
    except DataParseError as dfe:
        raise ValueError(str(dfe))
    tree_list = dataset.tree_lists[0]
    parent_id = '805080'
    branch_counter = 0
    tree_labels = set()
    for tree in tree_list:
        encode_splits(tree)
        tree_mask = tree.seed_node.edge.split_bitmask
        assert tree_mask is not None
        tree_tax = set(split_to_list(tree_mask))
        split_list = []
        for node in tree.leaf_iter():
            tree_labels.add(node.taxon.label)
    output.write('805080\t|\t\t|\tlife\t|\tno rank\t|\tncbi:1,gbif:0\t|\t\t|\t\n')
    for i in tree_labels:
        name, ottoid = i.split('@')
        output.write(ottoid+'\t|\t'+parent_id+'\t|\t'+name+'\t|\tspecies\t|\tncbi:1\n')

if __name__ == '__main__':
    try:
        filename = sys.argv[1]
    except:
        fo = sys.stdin
    else:
        fo = open(filename, 'rU')
    life_ott(fo)

