## Copyright (C) 2011, 2012 Robert Lipshitz and Sucharit Sarkar.
## Pure Python port (C) 2026.
## Licensed under GPLv3.

"""
convert.py -- Parse PD-presentation data from Knot Atlas .raw files
and write alldata.sage (knot definitions as CycList / ResolvedKnot objects).

Usage:  python convert.py
Input:  knots.raw, knots11.raw, links.raw
Output: alldata.sage
"""

from main import CycList, ResolvedKnot

outp = "alldata.sage"


def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


s = ''
for inpname in ['knots', 'knots11', 'links']:
    inp = inpname + '.raw'
    with open(inp, 'r') as reading:
        s = s + '\n' + reading.read()


s = (s.replace('<knot:', '')
      .replace('> <invariant:PD_Presentation> "X<sub>', ' ')
      .replace('</sub> X<sub>', ' ')
      .replace('</sub>" .', ''))
lines = s.split('\n')

while lines.count('') > 0:
    lines.remove('')


final = [None] * len(lines)
for i in range(len(lines)):
    final[i] = lines[i].split(' ', 1)
    final[i][0] = (final[i][0]
                   .replace('(', '_')
                   .replace(',', '_')
                   .replace(')', ''))
    if is_int(list(final[i][0])[0]):
        final[i][0] = 'R' + final[i][0]
    final[i][1] = final[i][1].split(' ')
    for j in range(len(final[i][1])):
        if final[i][1][j].find(',') == -1:
            final[i][1][j] = list(final[i][1][j])
        else:
            final[i][1][j] = final[i][1][j].split(',')
        for k in range(4):
            final[i][1][j][k] = int(final[i][1][j][k])


with open(outp, 'w') as writing:
    writing.write("all=[")
    for i in range(len(final)):
        writing.write("'" + final[i][0] + "',")
    writing.write("]\n")

    for i in range(len(final)):
        name = final[i][0]
        writing.write('\n#' + name + '\n')

        circles = []
        # flatten list of 4-element crossing lists into a flat list
        crossings = [x for sublist in final[i][1] for x in sublist]
        travelled = [False] * len(crossings)
        while False in travelled:
            temp = []
            complete = False
            start = travelled.index(False)
            now = start
            while not complete:
                travelled[now] = True
                if now % 2 == 0:
                    temp = temp + [-1 - int(now/4),]
                    now = now + 1
                else:
                    temp = temp + [1 + int(now/4),]
                    now = now - 1
                travelled[now] = True
                for j in range(len(crossings)):
                    if crossings[j] == crossings[now]:
                        if j == start:
                            complete = True
                        if not travelled[j]:
                            now = j
            circles = circles + [temp,]

        writing.write(name + '_def = [')
        for j in range(len(circles)):
            writing.write('CycList(' + repr(circles[j]) + '),')
        writing.write(']\n')

        [nplus, nminus] = [0, 0]
        for j in range(len(final[i][1])):
            [a, b] = [final[i][1][j][1], final[i][1][j][3]]
            if abs(a - b) == 1:
                if a > b:
                    nplus = nplus + 1
                if a < b:
                    nminus = nminus + 1
            else:
                if a > b:
                    nminus = nminus + 1
                if a < b:
                    nplus = nplus + 1
        writing.write(name + '_cross = ' + repr(tuple([nplus, nminus])) + '\n')
        writing.write(name + '_resolved_knot = ResolvedKnot(' + name + '_def)\n')
