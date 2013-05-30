import sys
class Taxon(object):
    def __init__(self,name,parent):
        self.name = name
        self.parent = parent
        self.children = []
        
    def write_as_newick(self, out):
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
            try:
                parent_id = int(row_list[1])
                parent = keyToObject.get(parent_id)
            except:
                parent = None
            taxon = Taxon(name, parent)
            keyToObject[t_id] = taxon
            if root is None:
                root = taxon
        
            if parent is not None:
               # print('parent of ' + t_id + ' is ' + parent.name)
                parent.children.append(taxon)
            else:
                pass
                #print('No parent for ' + name + 'looked up ' + str(parent_id))
                #print(str(keyToObject))
    return root
