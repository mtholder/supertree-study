class Taxon(object):
    def __init__(self,name,parent):
        self.name = name
        self.parent = parent
        self.children = []

keyToObject = {}
#root = none


def parse_ott(filename,splitchar ='\t|\t'):
    '''
    function to be used
    '''
    inp_file = open(filename, 'rU')
    inp = iter(inp_file)
    next(inp)
   # size = len(file.readlines())
    #print(size)
    for row in inp:
        row_list = row.split(splitchar)
        t_id = int(row_list[0])
        name = row_list[2]
        parent_id = int(row[1])
        parent = keyToObject.get(parent_id)
        taxon = Taxon(name, parent)
        keyToObject[t_id] = taxon
        if parent is not None:
            parent.children.append(taxon)
    inp_file.close()
