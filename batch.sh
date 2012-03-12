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


#!/bin/bash

echo 0 > stored.txt

OUT="Notdone"
while [[ "$OUT" == *Notdone* ]]
do
    OUT=`sage ./batch.sage`
done

rm stored.txt