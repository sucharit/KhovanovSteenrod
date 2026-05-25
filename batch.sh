#!/bin/bash

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

# Usage: ./batch.sh [--sage]
#   Default: use pure Python (python batch.py)
#   --sage:  use SageMath (sage batch.sage)
#
# Note: this script resets progress and starts from the beginning.
# To resume an interrupted run, omit this script and run the inner
# command in a loop directly:
#   while [[ "$(python batch.py)" == *Notdone* ]]; do :; done
#   while [[ "$(sage batch.sage)" == *Notdone* ]]; do :; done

if [[ "$1" == "--sage" ]]; then
    CMD="sage ./batch.sage"
else
    CMD="python batch.py"
fi

echo 0 > stored.txt

OUT="Notdone"
while [[ "$OUT" == *Notdone* ]]
do
    OUT=`$CMD`
done

rm stored.txt
