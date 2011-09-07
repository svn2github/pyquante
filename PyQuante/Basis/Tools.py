"""\
This program inputs a GAMESS basis, e.g. from the EMSL basis set order
form, and formats it into PyQuante format

Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

import re,pprint
from PyQuante.Element import name2no

basis_map = {
    '6-31g':'p631',
    '6-31g**':'p631ss',
    '6-31g(d,p)':'p631ss',
    '6-31g**++':'p631ppss',
    '6-31g++**':'p631ppss',
    '6-311g**':'p6311ss',
    '6-311g++(2d,2p)':'p6311pp_2d_2p',
    '6-311g++(3d,3p)':'p6311pp_3d_3p',
    '6-311g++(3df,3pd)':'p6311pp_3df_3pd',
    '3-21g':'p321',
    'sto3g':'sto3g',
    'sto-3g':'sto3g',
    'sto-6g':'sto6g',
    'lacvp':'lacvp',
    
    'ccpvdz':'ccpvdz',
    'cc-pvdz':'ccpvdz',
    'ccpvtz':'ccpvtz',
    'cc-pvtz':'ccpvtz',
    'ccpvqz':'ccpvqz',
    'cc-pvqz':'ccpvqz',
    'ccpv5z':'ccpv5z',
    'cc-pv5z':'ccpv5z',
    'ccpv6z':'ccpv6z',
    'cc-pv6z':'ccpv6z',

    'augccpvdz':'augccpvdz',
    'aug-cc-pvdz':'augccpvdz',
    'augccpvtz':'augccpvtz',
    'aug-cc-pvtz':'augccpvtz',
    'augccpvqz':'augccpvqz',
    'aug-cc-pvqz':'augccpvqz',
    'augccpv5z':'augccpv5z',
    'aug-cc-pv5z':'augccpv5z',
    'augccpv6z':'augccpv6z',
    'aug-cc-pv6z':'augccpv6z',
    
    'dzvp':'dzvp',
    }

def importname(modulename, name):
    """Import from a module whose name is determined at runtime.

    (Python Cookbook 2nd ed.)
    """
    module = __import__(modulename, globals(), locals(), [name])
    if not module:
        raise ImportError
    return getattr(module, name)

def get_basis_data(name):
    dc_name = name.lower()
    if dc_name not in basis_map:
        raise Exception("Can't import basis set %s %s" % (name,dc_name))
    return importname(basis_map[dc_name],"basis_data")

def split_comment(line):
    """Split a line into line,comment, where a comment is
    started by the ! character"""
    res = line.strip().split('!')
    if len(res) == 0: return "",""
    elif len(res) == 1: return res[0],""
    elif len(res) == 2: return res
    return res[0],''.join(res[1:])

def parse_gamess_basis(file,**kwargs):
    import re
    if type(file) == type(""): return parse_gamess_basis(open(file),**kwargs)
    maxatno = kwargs.get('maxatno',54)
    basis = {} # RPM changed basis sets from list to dictionary 2/22/2007
    atom_line = re.compile('[A-Za-z]{3,}')
    while 1:
        try:
            line = file.next()
        except:
            break
        #line,comment = split_comment(line)
        if line.startswith('!'):continue
        words = line.split()
        if not words: continue
        if atom_line.search(line):
            #el_string = line.strip()
            el_string = line.split()[0]
            atno = name2no[el_string]
            bfs = basis.setdefault(atno,[])
            #basis[atno] = bfs
        else:
            words = line.split()
            sym = words[0]
            nprim = int(words[1])
            if sym == "L":
                sprims = []
                pprims = []
                for i in xrange(nprim):
                    line = file.next()
                    words = line.split()
                    sprims.append((float(words[1]),float(words[2])))
                    pprims.append((float(words[1]),float(words[3])))
                bfs.append(("S",sprims))
                bfs.append(("P",pprims))
            else:
                prims = []
                for i in xrange(nprim):
                    line = file.next()
                    words = line.split()
                    prims.append((float(words[1]),float(words[2])))
                bfs.append((sym,prims))
    return basis

def main(**kwargs):
    fname = kwargs.get('fname','/home/rmuller/dzvp_basis.txt')
    oname = kwargs.get('oname','basis_dzvp.py')
    file = open(fname)
    basis = parse_gamess_basis(file,**kwargs)
    string = pprint.pformat(basis)
    file = open(oname,'w').write('basis_data = %s' % string)
    return

if __name__ == '__main__': main()
    
