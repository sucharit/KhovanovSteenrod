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

outp="allhom.sage"

def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


s=''
for inpname in ['knots', 'knots11', 'links']:
    inp=inpname+'hom.raw'
    reading = open(inp,'r')
    s=s+'\n'+reading.read()
    reading.close()


s=s.replace('<knot:','').replace('> <invariant:Integral_Khovanov_Homology> "{',' ').replace('\\n|}" .','')
lines=s.split('\n')

while lines.count('')>0:
    lines.remove('')


writing = open(outp,'w')

for i in range(len(lines)):
    [name,stuff]=lines[i].split(' ',1)

    name=name.replace('(','_').replace(',','_').replace(')','')
    if is_int(list(name)[0]):
        name='R'+name

    hom=dict()
    i_list=[]

    stuff=stuff.replace('/','').replace('<math>','').replace('\\\\','').replace('|','').replace('bgcolor=yellow','').replace('- align=center','XX').replace('{mathbb Z}','Z')
    stuff=stuff.replace('oplus',' ').replace('Z_2','F').replace('^','').replace('{','').replace('}','').split('\\n')
    while stuff.count('XX')>0:
        stuff.remove('XX')
    stuff=stuff[2:]

    
    for x in stuff:
        if list(x)[:2]==['i','=']:
            i_list=i_list+[int(x.lstrip('i=')),]

    stuff=stuff[len(i_list):]
    count=0
    for x in stuff:

        if count == 0:
            r = int(x.lstrip('r='))# r is hom grading
        else:
            i=i_list[count-1]
            q = 2*r + i # q is quantum (internal) grading

            x = x.split(' ')
            for y in x:
                if y !='':
                    typ=list(y)[0] # Z or F=Z/2?
                    y=y.lstrip(typ)
                    if y=='':
                        rank=1
                    else:
                        rank=int(y)

                    if not (r,q) in hom.keys():
                        hom[(r,q)]=0
                    hom[(r,q)]=hom[(r,q)]+rank

                    if typ=='F':
                        if not (r-1,q) in hom.keys():
                            hom[(r-1,q)]=0
                        hom[(r-1,q)]=hom[(r-1,q)]+rank
                
            
        count = (count+1)%(1+len(i_list))
        
    total_rank=[]
    for x in hom.keys():
        total_rank.append((x,hom[x]))
    total_rank.sort(key=lambda x: (x[0][1],x[0][0]) )

    writing.write('\n'+name+'_hom='+repr(total_rank)+'\n')

writing.close()

