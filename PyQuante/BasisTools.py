#!/usr/bin/env python
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
    maxatno = kwargs.get('maxatno',54)
    basis = [None]*(maxatno+1)
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
            el_string = line.strip()
            atno = name2no[el_string]
            bfs = []
            basis[atno] = bfs
        else:
            words = line.split()
            sym = words[0]
            nprim = int(words[1])
            prims = []
            for i in range(nprim):
                line = file.next()
                words = line.split()
                prims.append((float(words[1]),float(words[2])))
            bfs.append((sym,prims))
    return basis

def main(**opts):
    fname = opts.get('fname','/home/rmuller/dzvp_basis.txt')
    oname = opts.get('oname','basis_dzvp.py')
    file = open(fname)
    basis = parse_gamess_basis(file,**opts)
    string = pprint.pformat(basis)
    file = open(oname,'w').write('basis_data = %s' % string)
    return

if __name__ == '__main__': main()
    
