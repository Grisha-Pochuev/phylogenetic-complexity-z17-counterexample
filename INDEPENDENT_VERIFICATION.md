# Independent verification protocol

This file gives a minimal path for checking the counterexample without relying on the saved output files.

## 1. Check the primary verifier

Run from the repository root:

```bash
python verification/verify_phylogenetic_Z17.py
```

The script uses only the Python standard library. It reconstructs all allowed row types by checking all \(17^3\) triples, verifies zero row sums and column multiplicities, builds the incidence matrix with exact rational arithmetic, verifies rank 21 and the determinant-17 minor certificate, derives the one-dimensional kernel, and independently enumerates the nonnegative compatible count vectors.

Expected final conclusions include:

```text
rank 21
nullity 1
21 x 21 minor determinant 17
exactly 2 compatible tables
minimum nontrivial move size 18
```

## 2. Check with a different algorithm

Run:

```bash
python verification/audit_Z17_independent.py
```

This script does not use matrix rank, the determinant certificate, or the one-dimensional kernel conclusion. It fixes the first-column entries in sorted order and counts completions from the remaining multiplicities of the second and third columns by dynamic programming.

Expected JSON values:

```json
{
  "labelled_fillings_total": 96,
  "labelled_fillings_A": 24,
  "labelled_fillings_B": 72,
  "dynamic_programming_states": 981,
  "unordered_compatible_tables": 2,
  "minimum_nontrivial_move_size": 18
}
```

The identity \(96=24+72\) leaves no room for a third compatible unordered table.

## 3. Optional symbolic audit

Install the exact version used for the stored symbolic run:

```bash
python -m pip install sympy==1.14.0
python verification/audit_Z17_symbolic.py
```

Expected symbolic invariants:

- incidence matrix shape: \(23\times22\);
- rank over \(\mathbb Q\): 21;
- the displayed kernel basis is \(B-A\);
- the specified maximal minor has determinant 17.

## 4. Read the proof independently of the code

The short English manuscript is in `paper/`. The expanded Russian proof is in `proof/`. The mathematical core is elementary: once the compatible fiber is proved to contain only \(A\) and \(B\), disjoint row supports imply that a move changing at most 17 of 18 rows cannot leave \(A\) nontrivially.

## 5. Manual-only GitHub Actions

`.github/workflows/verify.yml` has only a `workflow_dispatch` trigger. It will not run when files are pushed. A reviewer may run it manually in a fork, which uses the fork's Actions allocation rather than relying on the author's currently occupied runners.

## 6. What would falsify the package

Any one of the following would invalidate the claimed certificate:

- a row in either table fails to sum to 0 modulo 17;
- the corresponding column multisets differ;
- an allowed row type inside the three supports is missing from the list of 22;
- a third nonnegative integer count vector has the same marginals;
- the two tables share a row type;
- the incidence matrix has rank different from 21;
- the stated minor determinant is not 17.

The two core programs test these points by substantially different routes.
