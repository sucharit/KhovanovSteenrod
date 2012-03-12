## Copyright (C) 2011, 2012 Robert Lipshitz and Sucharit Sarkar.
## Contact: lipshitz@math.columbia.edu, sucharit@math.columbia.edu

## This file is part of KhovanovSteenrod.

## KhovanovSteenrod is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.

## KhovanovSteenrod is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with KhovanovSteenrod; see COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

outp="alldata.sage"

def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


s=''
for inpname in ['knots', 'knots11', 'links']:
    inp=inpname+'.raw'
    reading = open(inp,'r')
    s=s+'\n'+reading.read()
    reading.close()


s=s.replace('<knot:','').replace('> <invariant:PD_Presentation> "X<sub>',' ').replace('</sub> X<sub>',' ').replace('</sub>" .','')
lines=s.split('\n')

while lines.count('')>0:
    lines.remove('')


final=[None]*(len(lines))
for i in range(len(lines)):
    final[i]=lines[i].split(' ',1)
    final[i][0]=final[i][0].replace('(','_').replace(',','_').replace(')','')
    if is_int(list(final[i][0])[0]):
        final[i][0]='R'+final[i][0]
    final[i][1]=final[i][1].split(' ')
    for j in range(len(final[i][1])):
        if final[i][1][j].find(',') == -1:
            final[i][1][j]=list(final[i][1][j])
        else:
            final[i][1][j]=final[i][1][j].split(',')
        for k in range(4):
            final[i][1][j][k]=int(final[i][1][j][k])



writing = open(outp,'w')
writing.write("all=[")
for i in range(len(final)):
    writing.write('\''+final[i][0]+'\',')
writing.write("]\n")


for i in range(len(final)):
    name=final[i][0]
    writing.write('\n#'+name+'\n')

    circles=[]
    crossings=flatten(final[i][1])
    travelled=[False,]*len(crossings)
    while False in travelled:
        temp=[]
        complete=False
        start=travelled.index(False)
        now = start
        while not complete:
            travelled[now]=True
            if now%2 == 0:
                temp = temp + [-1-int(now/4),]
                now=now+1;
            else:
                temp = temp + [1+int(now/4),]
                now=now-1;
            travelled[now]=True
            for j in range(len(crossings)):
                if crossings[j] == crossings[now]:
                    if j == start:
                        complete = True
                    if not travelled[j]:
                        now = j
        circles=circles+[temp,]
    writing.write(name+'_def = [')
    for j in range(len(circles)):
        writing.write('CycList('+repr(circles[j])+'),')
    writing.write(']\n')

    [nplus,nminus]=[0,0]
    for j in range(len(final[i][1])):
        [a,b]=[final[i][1][j][1],final[i][1][j][3]]
        if abs(a-b)==1:
            if a>b:
                nplus=nplus+1
            if a<b:
                nminus=nminus+1
        else:
            if a>b:
                nminus=nminus+1
            if a<b:
                nplus=nplus+1
    writing.write(name+'_cross = '+repr(tuple([nplus,nminus]))+'\n')

    writing.write(name+'_resolved_knot = ResolvedKnot('+name+'_def)\n')


writing.close()

