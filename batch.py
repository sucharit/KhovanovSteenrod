## Copyright (C) 2011, 2012 Robert Lipshitz and Sucharit Sarkar.
## Pure Python port (C) 2026.
## Licensed under GPLv3.

"""
batch.py -- Compute Sq^1 and Sq^2 for all knots/links with checkpointing.

Usage:  python batch.py          (or:  nohup python batch.py &)
Input:  alldata.sage, allhom.sage
Output: result.sage  (Steenrod square matrices, appended incrementally)
        stored.txt   (checkpoint: index of next knot to process)

The script can be interrupted and restarted; it resumes from where it left off.
"""

import sys
from main import CycList, ResolvedKnot, KnotKh

outp   = 'result.sage'
stored = 'stored.txt'

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

_ns = {'CycList': CycList, 'ResolvedKnot': ResolvedKnot}
with open('alldata.sage') as _f:
    exec(_f.read(), _ns)

with open('allhom.sage') as _f:
    exec(_f.read(), _ns)

univ_list   = _ns['all']
univ_number = len(univ_list)

# ---------------------------------------------------------------------------
# Read checkpoint
# ---------------------------------------------------------------------------

start = 0
try:
    with open(stored, 'r') as _data:
        s = _data.read().strip()
        if s:
            start = int(s)
except FileNotFoundError:
    pass

# ---------------------------------------------------------------------------
# Open output file
# ---------------------------------------------------------------------------

if start == 0:
    writing = open(outp, 'w')
    writing.write('all=' + repr(univ_list) + '\n')
    writing.flush()
else:
    writing = open(outp, 'a')

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

written_once = False

for i in range(start, univ_number):
    K   = univ_list[i]
    hom = dict(_ns[K + '_hom'])

    writing.write('\n\n' + K + '_hom=' + repr(hom))
    writing.flush()

    sq1 = dict()
    sq2 = dict()
    knot_defined = False

    for (t, q) in hom.keys():
        # Only compute Sq1/Sq2 in bigradings where the result is potentially nontrivial
        compute_sq2 = (t+2, q) in hom
        compute_sq1_only = ((t+1, q) in hom and (t-1, q) in hom)

        if compute_sq2 or compute_sq1_only:
            written_once = True
            if not knot_defined:
                resolved_knot = _ns[K + '_resolved_knot']
                cross         = _ns[K + '_cross']
                Central_Knot  = KnotKh(resolved_knot, cross, print_progress=False)
                knot_defined  = True

            if compute_sq2:
                sq1[(t, q)] = Central_Knot.sq1((t, q), print_progress=False).rows()
                sq2[(t, q)] = Central_Knot.sq2((t, q), print_progress=False).rows()
            elif compute_sq1_only:
                sq1[(t, q)] = Central_Knot.sq1((t, q), print_progress=False).rows()

    writing.write('\n' + K + '_sq1=' + repr(sq1))
    writing.write('\n' + K + '_sq2=' + repr(sq2))
    writing.flush()

    if written_once:
        with open(stored, 'w') as storing:
            storing.write(repr(i + 1))
        break

if written_once:
    print("Notdone")
else:
    print("Done")

writing.close()
