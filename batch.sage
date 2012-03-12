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

load main.sage
load alldata.sage
load allhom.sage

univ_list = all
univ_number = len(univ_list)

outp = 'result.sage'
stored = 'stored.txt'


written_once=False


data=open(stored,'r')
s=data.read()
if s!='':
    n=int(s)
    start=n



if start == 0:
    writing = open(outp,'w')
    writing.write('all='+repr(univ_list)+'\n')
    writing.flush()

else:
    writing = open(outp,'a')



for i in range(start,univ_number):

    K = univ_list[i]
    hom = dict(eval(K+'_hom'))

    writing.write('\n\n'+K+'_hom='+repr(hom))
    writing.flush()

    sq1=dict()
    sq2=dict()
    knot_defined=False
    
    for (t,q) in hom.keys():#Only computing Sq1 in wide rows.
        if (t+2,q) in hom.keys():
            written_once=True
            if not knot_defined:
                compute_knot(K)
                knot_defined=True
            sq1[(t,q)]=Central_Knot.sq1((t,q),print_progress=False).rows()
            sq2[(t,q)]=Central_Knot.sq2((t,q),print_progress=False).rows()
        elif (t+1,q) in hom.keys() and (t-1,q) in hom.keys():
            written_once=True
            if not knot_defined:
                compute_knot(K)
                knot_defined=True
            sq1[(t,q)]=Central_Knot.sq1((t,q),print_progress=False).rows()

    writing.write('\n'+K+'_sq1='+repr(sq1))
    writing.write('\n'+K+'_sq2='+repr(sq2))
    
    writing.flush()
    if written_once:
        storing = open(stored,'w')
        storing.write(repr(i+1))
        break

if written_once:
    print "Notdone"
else:
    print "Done"

writing.close()
