#!/usr/bin/env python
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
        '''
        This function writes the output in form of Newick.
        '''
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
    #taxon_list = []
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
            mod_name = mod_name.replace("/","''/''")
            #print mod_name
            #if name is not "Microcebus sp. d'Ambre":
            #name.replace("'","''")
            try:
                parent_id = int(row_list[1])
                parent = keyToObject.get(parent_id)
            except:
                parent = None      
            taxon = Taxon(mod_name, parent)
            keyToObject[t_id] = taxon
            #taxon_list.append(taxon) #temperary
            if root is None:
                root = taxon
            if parent is not None:
                #print('parent of ' + t_id + ' is ' + parent.name)
                parent.children.append(taxon)
            else:
                pass
                #print('No parent for ' + name + 'looked up ' + str(parent_id))
                #print(str(keyToObject))   
    # for t in taxon_list:
    #     if t.parent is None:
    #         print "Parent list taxon:", t.name       
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
    #primate_id_set = set()
    virdi_id_set = set()
    with open(filename, 'rU') as inp_file:
        inp = iter(inp_file)
        #next(inp)
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
                #if name == 'Primates' or parent_id in primate_id_set:
                if name == 'Chordata' or parent_id in virdi_id_set:
                    sys.stdout.write(row)
                    #primate_id_set.add(t_id)
                    virdi_id_set.add(t_id)


if __name__ == '__main__':
   r = parse_ott(sys.argv[1])
   r.write_as_newick(sys.stdout)
   '''
   print '\n'
   e = parse_ott2(sys.argv[1])
   e.write_as_newick(sys.stdout)
   '''


