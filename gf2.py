## Copyright (C) 2011, 2012 Robert Lipshitz and Sucharit Sarkar.
## Pure Python port (C) 2026.
## Licensed under GPLv3.

"""GF(2) linear algebra using numpy.

Provides GF2Matrix (matrix over GF(2)) and helper functions for
Gaussian elimination, kernel, row space, rank, and linear system
solving over the field with 2 elements.
"""

import numpy as np


# ---------------------------------------------------------------------------
# Low-level algorithms (operate on numpy uint8 2D arrays)
# ---------------------------------------------------------------------------

def _rref(M):
    """Reduced row echelon form of M over GF(2).

    Returns (rref_M, pivot_cols, rank).
    """
    M = M.copy().astype(np.uint8)
    nrows, ncols = M.shape
    pivot_row = 0
    pivot_cols = []
    for col in range(ncols):
        # Find a pivot in this column at or below pivot_row
        found = -1
        for row in range(pivot_row, nrows):
            if M[row, col]:
                found = row
                break
        if found == -1:
            continue
        # Swap rows
        if found != pivot_row:
            M[[pivot_row, found]] = M[[found, pivot_row]]
        # Eliminate all other rows (full RREF, not just REF)
        for row in range(nrows):
            if row != pivot_row and M[row, col]:
                M[row] ^= M[pivot_row]
        pivot_cols.append(col)
        pivot_row += 1
    return M, pivot_cols, pivot_row


def _rank(M):
    """Rank of M over GF(2)."""
    if M.size == 0:
        return 0
    _, _, r = _rref(M)
    return r


def _row_space(M):
    """Basis of row space of M over GF(2).

    Returns list of numpy uint8 1-D arrays (the nonzero rows of RREF(M)).
    """
    if M.size == 0:
        return []
    rref, _, rank = _rref(M)
    return [rref[i].copy() for i in range(rank)]


def _right_null_space(M):
    """Right null space of M: {x : M @ x = 0} over GF(2).

    Returns list of numpy uint8 1-D arrays.
    """
    if M.size == 0:
        return []
    nrows, ncols = M.shape
    if ncols == 0:
        return []
    if nrows == 0:
        return [np.eye(ncols, dtype=np.uint8)[k] for k in range(ncols)]

    rref, pivot_cols, rank = _rref(M)
    pivot_set = set(pivot_cols)
    free_cols = [c for c in range(ncols) if c not in pivot_set]

    null_basis = []
    for fc in free_cols:
        x = np.zeros(ncols, dtype=np.uint8)
        x[fc] = 1
        for i, pc in enumerate(pivot_cols):
            x[pc] = rref[i, fc]   # in GF(2): -entry = entry
        null_basis.append(x)
    return null_basis


def _left_null_space(M):
    """Left null space of M: {v : v @ M = 0} over GF(2).

    Returns list of numpy uint8 1-D arrays.
    """
    if M.size == 0:
        return []
    return _right_null_space(M.T.copy())


def _pivot_rows(M):
    """Indices of rows of M that contribute to its row space (over GF(2)).

    Processes rows in order; a row is a 'pivot row' if it is not already
    in the span of the rows before it.
    """
    if M.size == 0:
        return []
    nrows, ncols = M.shape
    pivot_col_map = {}   # pivot_column -> current basis row (numpy array)
    pivot_row_indices = []

    for orig_idx in range(nrows):
        row = M[orig_idx].astype(np.uint8).copy()
        # Reduce row against current basis vectors
        for col in range(ncols):
            if row[col] and col in pivot_col_map:
                row ^= pivot_col_map[col]
        # If the reduced row is nonzero, this row is a new pivot row
        nz = np.nonzero(row)[0]
        if len(nz):
            lead = int(nz[0])
            pivot_col_map[lead] = row
            pivot_row_indices.append(orig_idx)

    return pivot_row_indices


def _solve_left(A, b):
    """Find x (1-D numpy uint8 array) such that x @ A = b over GF(2).

    Equivalent to solving A^T @ x^T = b^T.
    Returns one solution (free variables set to 0).
    Raises ValueError if the system is inconsistent.
    """
    nrows, ncols = A.shape   # x has length nrows, b has length ncols
    if nrows == 0:
        return np.zeros(0, dtype=np.uint8)

    At = A.T.astype(np.uint8).copy()   # shape: ncols x nrows
    b_arr = b.astype(np.uint8).copy()  # length: ncols

    # Augmented system [At | b_arr], shape ncols x (nrows + 1)
    aug = np.hstack([At, b_arr.reshape(-1, 1)])

    n_eqs = aug.shape[0]   # = ncols
    n_vars = nrows

    pivot_col_to_row = {}   # variable index -> pivot row index
    pivot_row_idx = 0
    for col in range(n_vars):
        found = -1
        for row in range(pivot_row_idx, n_eqs):
            if aug[row, col]:
                found = row
                break
        if found == -1:
            continue
        if found != pivot_row_idx:
            aug[[pivot_row_idx, found]] = aug[[found, pivot_row_idx]]
        for row in range(n_eqs):
            if row != pivot_row_idx and aug[row, col]:
                aug[row] ^= aug[pivot_row_idx]
        pivot_col_to_row[col] = pivot_row_idx
        pivot_row_idx += 1

    # Check consistency
    for row in range(pivot_row_idx, n_eqs):
        if aug[row, n_vars]:
            raise ValueError("Inconsistent system in _solve_left")

    # Extract solution (free variables default to 0)
    x = np.zeros(n_vars, dtype=np.uint8)
    for col, pidx in pivot_col_to_row.items():
        x[col] = aug[pidx, n_vars]
    return x


# ---------------------------------------------------------------------------
# GF2Matrix class
# ---------------------------------------------------------------------------

class GF2Matrix:
    """Matrix over GF(2), backed by a numpy uint8 2-D array.

    Creation:
        GF2Matrix(nrows, ncols)                  -- zero matrix
        GF2Matrix(nrows, ncols, flat_data)       -- from flat list
        GF2Matrix(nrows, ncols, 2d_array)        -- from 2-D array
        GF2Matrix.from_rows(list_of_rows)        -- from list of row vectors

    Indexing:
        M[i, j]   -> int (0 or 1)
        M[i, j] = val   -> set entry (automatically reduced mod 2)
        M[i]      -> numpy uint8 1-D array (copy of i-th row)
        M[i] = v  -> set i-th row

    Operations:
        M.rows()           -> list of tuples (for serialisation)
        M.transpose()      -> new GF2Matrix
        M.rank()           -> int
        M.kernel()         -> GF2RowSpace (left null space)
        M.row_space()      -> GF2RowSpace
        M.pivot_rows()     -> list of int
        M.solve_left(b)    -> numpy uint8 1-D array x with x @ M = b
        A @ B              -> GF2Matrix or numpy array
        v @ M  (numpy v)   -> numpy uint8 1-D array
    """

    def __init__(self, nrows, ncols, data=None):
        self.nrows = nrows
        self.ncols = ncols
        if data is None:
            self._m = np.zeros((nrows, ncols), dtype=np.uint8)
        else:
            arr = np.asarray(data, dtype=np.int64)
            if arr.ndim == 1:
                self._m = (arr % 2).astype(np.uint8).reshape(nrows, ncols)
            else:
                self._m = (arr % 2).astype(np.uint8)

    @classmethod
    def from_rows(cls, rows):
        """Create matrix from a list of row vectors (lists/tuples/arrays/ints).

        Handles Sage's serialisation format where a 1-element GF(2) vector is
        stored as a bare integer (0 or 1) rather than as a 1-tuple.
        """
        if not rows:
            return cls(0, 0)
        # Normalise each row: bare int → 1-element list
        normalised = []
        for r in rows:
            if isinstance(r, (int, np.integer)):
                normalised.append(np.array([int(r)], dtype=np.uint8) % 2)
            else:
                normalised.append(np.asarray(r, dtype=np.uint8) % 2)
        ncols = len(normalised[0])
        nrows = len(normalised)
        obj = cls(nrows, ncols)
        for i, r in enumerate(normalised):
            obj._m[i, :len(r)] = r
        return obj

    def __getitem__(self, key):
        if isinstance(key, tuple):
            i, j = key
            return int(self._m[i, j])
        else:
            return self._m[key].copy()   # numpy array (i-th row)

    def __setitem__(self, key, val):
        if isinstance(key, tuple):
            i, j = key
            self._m[i, j] = int(val) % 2
        else:
            self._m[key] = np.asarray(val, dtype=np.uint8) % 2

    def rows(self):
        """Return list of rows as tuples (each entry is 0 or 1)."""
        return [tuple(int(x) for x in self._m[i]) for i in range(self.nrows)]

    def transpose(self):
        """Return the transpose as a new GF2Matrix."""
        obj = GF2Matrix(self.ncols, self.nrows)
        obj._m = self._m.T.copy()
        return obj

    def rank(self):
        """Rank of this matrix over GF(2)."""
        return _rank(self._m)

    def kernel(self):
        """Left null space as a GF2RowSpace: {v : v @ self = 0}."""
        return GF2RowSpace(_left_null_space(self._m), self.nrows)

    def kernel_basis(self):
        """Left null space basis as a list of numpy uint8 arrays."""
        return _left_null_space(self._m)

    def row_space(self):
        """Row space as a GF2RowSpace."""
        return GF2RowSpace(_row_space(self._m), self.ncols)

    def row_space_basis(self):
        """Row space basis as a list of numpy uint8 arrays."""
        return _row_space(self._m)

    def pivot_rows(self):
        """Indices of rows that contribute to the row space."""
        return _pivot_rows(self._m)

    def solve_left(self, b):
        """Find x such that x @ self = b over GF(2). Returns numpy uint8 array."""
        b_arr = np.asarray(b, dtype=np.uint8) % 2
        return _solve_left(self._m, b_arr)

    def __matmul__(self, other):
        """self @ other."""
        if isinstance(other, GF2Matrix):
            result = self._m.astype(np.int32) @ other._m.astype(np.int32) % 2
            return GF2Matrix(self.nrows, other.ncols, result)
        elif isinstance(other, np.ndarray):
            return (self._m.astype(np.int32) @ other.astype(np.int32) % 2).astype(np.uint8)

    def __rmatmul__(self, other):
        """other @ self  (when other is a numpy array)."""
        if isinstance(other, np.ndarray):
            return (other.astype(np.int32) @ self._m.astype(np.int32) % 2).astype(np.uint8)

    def __repr__(self):
        return f"GF2Matrix({self.nrows}x{self.ncols})"


# ---------------------------------------------------------------------------
# GF2RowSpace class  (returned by GF2Matrix.row_space() and .kernel())
# ---------------------------------------------------------------------------

class GF2RowSpace:
    """Represents a subspace of GF(2)^n given by a list of basis vectors."""

    def __init__(self, basis_rows, ncols=None):
        """
        basis_rows: list of numpy uint8 1-D arrays.
        ncols:      ambient dimension (inferred from basis if not given).
        """
        self._basis = list(basis_rows)
        if ncols is not None:
            self._ncols = ncols
        elif basis_rows:
            self._ncols = len(basis_rows[0])
        else:
            self._ncols = 0

    def basis(self):
        """Return basis as a list of numpy uint8 arrays."""
        return list(self._basis)

    def rank(self):
        """Dimension of this subspace."""
        return len(self._basis)

    def intersection(self, other):
        """Return GF2RowSpace representing the intersection with another subspace.

        Uses the dimension formula:
            dim(A ∩ B) = dim(A) + dim(B) - dim(A + B)
        and explicitly computes a basis via the left null space of [A; B].
        """
        if not self._basis or not other._basis:
            return GF2RowSpace([], self._ncols)

        A = np.array(self._basis, dtype=np.uint8)
        B = np.array(other._basis, dtype=np.uint8)
        rA, rB = len(self._basis), len(other._basis)

        # Find [x; y] with x @ A = y @ B  (equivalently [x,y] @ [A;B] = 0 in GF(2))
        AB = np.vstack([A, B])   # (rA+rB) x ncols
        null_vecs = _left_null_space(AB)  # each vec has length rA+rB

        inter_vecs = []
        for v in null_vecs:
            x = v[:rA]
            if not np.any(x):
                continue
            w = (x.astype(np.int32) @ A.astype(np.int32) % 2).astype(np.uint8)
            if np.any(w):
                inter_vecs.append(w)

        if not inter_vecs:
            return GF2RowSpace([], self._ncols)

        # Reduce to a basis
        basis = _row_space(np.array(inter_vecs, dtype=np.uint8))
        return GF2RowSpace(basis, self._ncols)
