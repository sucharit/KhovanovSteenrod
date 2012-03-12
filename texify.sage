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

load result.sage

writing=open('out.tex','w')

for K in all:
    is_trivial=True
    
    hom=eval(K+'_hom')
    sq1=eval(K+'_sq1')
    sq2=eval(K+'_sq2')
    st = dict()

    for (t,q) in sq2.keys():
        whole=matrix(GF(2),hom[(t,q)],hom[(t+2,q)],sq2[(t,q)])

        if (t+1,q) in hom.keys():#Minor hack
            tail=matrix(GF(2),hom[(t,q)],hom[(t+1,q)],sq1[(t,q)])
            head=matrix(GF(2),hom[(t+1,q)],hom[(t+2,q)],sq1[(t+1,q)])
        else:
            tail=matrix(GF(2),hom[(t,q)],0)
            head=matrix(GF(2),0,hom[(t+2,q)])

        a = whole.rank()
        if a>0:
            is_trivial=False

            st[(t,q)]=[0,0,0,0]
            st[(t,q)][0]=a

            ker_basis=tail.kernel().basis()
            restricted=matrix(GF(2),[v*whole for v in ker_basis])

            st[(t,q)][1]=restricted.rank()
            
            im=head.row_space()

            if im.rank()>0:
                st[(t,q)][2]=whole.row_space().intersection(im).rank()
                if st[(t,q)][1]>0:
                    st[(t,q)][3]=restricted.row_space().intersection(im).rank()

            [r1,r2,r3,r4]=st[(t,q)]
            st[(t,q)]=[r2-r4,r1-r2-r3+r4,r4,r3-r4]


    if not is_trivial:
        K_first=K[0]
        if K_first == 'R':
            K=K.lstrip('R').replace('_','_{')+'}'
        if K_first == 'K':
            K=K.replace('K','\\text{K}')
        if K_first == 'L':
            K=K.replace('L','\\text{L}')
        writing.write('$'+K+'$ & ')

        first_item=True
        for x in st.keys():
            if not first_item:
                writing.write(', ')
            else:
                first_item=False
            writing.write('$'+repr(x)+'\\mapsto '+repr(tuple(st[x]))+'$')
        writing.write('\\\\\n')
        writing.flush()
            
    

writing.close()
