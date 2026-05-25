## Copyright (C) 2011, 2012 Robert Lipshitz and Sucharit Sarkar.
## Pure Python port (C) 2026.
## Licensed under GPLv3.

"""
extract.py -- Parse Knot Atlas integral Khovanov homology data from .raw files
and write allhom.sage.

Usage:  python extract.py
Input:  knotshom.raw, knots11hom.raw, linkshom.raw
Output: allhom.sage
"""

outp = "allhom.sage"


def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


s = ''
for inpname in ['knots', 'knots11', 'links']:
    inp = inpname + 'hom.raw'
    with open(inp, 'r') as reading:
        s = s + '\n' + reading.read()


s = (s.replace('<knot:', '')
      .replace('> <invariant:Integral_Khovanov_Homology> "{', ' ')
      .replace('\\n|}" .', ''))
lines = s.split('\n')

while lines.count('') > 0:
    lines.remove('')


with open(outp, 'w') as writing:
    for i in range(len(lines)):
        [name, stuff] = lines[i].split(' ', 1)

        name = name.replace('(', '_').replace(',', '_').replace(')', '')
        if is_int(list(name)[0]):
            name = 'R' + name

        hom = dict()
        i_list = []

        stuff = (stuff.replace('/', '')
                      .replace('<math>', '')
                      .replace('\\\\', '')
                      .replace('|', '')
                      .replace('bgcolor=yellow', '')
                      .replace('- align=center', 'XX')
                      .replace('{mathbb Z}', 'Z'))
        stuff = (stuff.replace('oplus', ' ')
                      .replace('Z_2', 'F')
                      .replace('^', '')
                      .replace('{', '')
                      .replace('}', '')
                      .split('\\n'))
        while stuff.count('XX') > 0:
            stuff.remove('XX')
        stuff = stuff[2:]

        for x in stuff:
            if list(x)[:2] == ['i', '=']:
                i_list = i_list + [int(x.lstrip('i=')),]

        stuff = stuff[len(i_list):]
        count = 0
        for x in stuff:
            if count == 0:
                r = int(x.lstrip('r='))
            else:
                i_val = i_list[count-1]
                q = 2*r + i_val

                x = x.split(' ')
                for y in x:
                    if y != '':
                        typ = list(y)[0]   # 'Z' or 'F' (=Z/2)
                        y = y.lstrip(typ)
                        if y == '':
                            rank = 1
                        else:
                            rank = int(y)

                        if (r, q) not in hom:
                            hom[(r, q)] = 0
                        hom[(r, q)] += rank

                        if typ == 'F':
                            if (r-1, q) not in hom:
                                hom[(r-1, q)] = 0
                            hom[(r-1, q)] += rank

            count = (count + 1) % (1 + len(i_list))

        total_rank = list(hom.items())
        total_rank.sort(key=lambda x: (x[0][1], x[0][0]))

        writing.write('\n' + name + '_hom=' + repr(total_rank) + '\n')
