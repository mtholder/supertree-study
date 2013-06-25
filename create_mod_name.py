import sys
import string
if __name__ == '__main__':
    '''
    simple method to take the ottolid 
    and '@' symbol off of the taxa name 
    for the MRP algorithm to compute
    '''
    output = sys.stdout
    fo = open(sys.argv[1], "rU")
    for i in fo:
        taxa = i.translate(None, string.digits)
        taxa = taxa.translate(None, '@')
        print taxa

