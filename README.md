# KhovanovSteenrod: computing the first two Steenrod squares on Khovanov homology

This is the README file for the KhovanovSteenrod program, code for computing the first two Steenrod squares on Khovanov homology, induced by the Khovanov stable homotopy type.

## Usage

Here is a step by step guideline on how to recompute the data. There are two options: using Sage or using pure Python. The functionality is essentially identical.

This process will regenerate the files `alldata.sage`, `allhom.sage`, and `result.sage`.

- If using Sage, install [Sage](http://www.sagemath.org/).

- If instead using Python, install [numpy](https://numpy.org/):
  ```bash
  pip install numpy
  ```

- Open a terminal and go to the directory where you have downloaded all
the files of KhovanovSteenrod.  Henceforth, we are assuming that you
are using a Bash shell (the other shells will be similar), and Bash is
configured so that the command `sage` runs the program Sage.

- Get the data tables from Knot Atlas.
  ```bash
  wget http://katlas.org/Data/Rolfsen.rdf.gz
  wget http://katlas.org/Data/Knots11.rdf.gz
  wget http://katlas.org/Data/Links.rdf.gz
  ```
  On macOS, `wget` is not installed by default; replace each `wget <url>` with `curl -O <url>`.

- Extract the PD-presentations.
  ```bash
  zgrep PD_Presentation Rolfsen.rdf.gz > knots.raw
  zgrep PD_Presentation Knots11.rdf.gz > knots11.raw
  zgrep PD_Presentation Links.rdf.gz > links.raw
  ```

- Remove the unknot and the Hopf link. (Due to an unfortunate premature
  optimization, `main.sage`/`main.py` does not work if the link
  diagram has an unknot or Hopf link component.)
  ```bash
  sed -i 1d knots.raw    # Linux
  sed -i 1d links.raw    # Linux
  ```
  On macOS, use `sed -i '' 1d` instead of `sed -i 1d`.

- Convert the PD-presentations to our format. (This generates the file
  `alldata.sage`.)
  ```bash
  sage convert.sage   # Sage
  python convert.py   # pure Python
  ```
  (You only need to run one of these two commands; and similarly in later steps.)

- Compute Kh (over $F_2$) for all prime links up to 11 crossings. (This
  generates the file `allhom.sage`.) This can be done in two ways.

    - Compute the homology directly; this takes a long time.
      ```bash
      sage computehom.sage   # Sage
      python computehom.py   # pure Python
      ```

    - Extract from the Knot Atlas data.
      ```bash
      zgrep Integral_Khovanov Rolfsen.rdf.gz > knotshom.raw
      zgrep Integral_Khovanov Knots11.rdf.gz > knots11hom.raw
      zgrep Integral_Khovanov Links.rdf.gz > linkshom.raw
      sed -i 1d linkshom.raw   # Linux; macOS: sed -i '' 1d linkshom.raw
      sage extract.sage   # Sage
      python extract.py   # pure Python
      ```

- Compute the action of $\mathrm{Sq}^1$ and $\mathrm{Sq}^2$ on all quantum gradings where the
  width is at least 3. (This generates the file `result.sage`.)
  This process might take forever, so we recommend using `nohup`.
  ```bash
  nohup ./batch.sh &           # pure Python (default)
  nohup ./batch.sh --sage &    # SageMath
  ```
  **Note:** `batch.sh` always restarts from the beginning. To resume an
  interrupted run without losing progress, run the inner loop directly:
  ```bash
  while [[ "$(python batch.py)" == *Notdone* ]]; do :; done       # pure Python
  while [[ "$(sage ./batch.sage)" == *Notdone* ]]; do :; done     # SageMath
  ```

- Generate the LaTeX snippet from `result.sage`.
  ```bash
  sage texify.sage   # Sage
  python texify.py   # pure Python
  ```

## Credits, copyright, and license

Copyright (c) 2011 - 2026 by Robert Lipshitz and Sucharit Sarkar.
Contact: lipshitz@uoregon.edu, sucharit@math.ucla.edu

Python port created by Anthropic's Claude chatbot (Sonnet 4.6), with mild assistance by Robert Lipshitz.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
