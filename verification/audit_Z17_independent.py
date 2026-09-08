#!/usr/bin/env python3
"""Independent exact audit of the Z/17Z counterexample.

Counts labelled second/third-column fillings with the first column fixed.
Uses dynamic programming on remaining multiplicities, rather than the
linear-algebra / grouped-multiplicity methods in verify_phylogenetic_Z17.py.
Requires Python 3.9+ and its standard library only.
"""
from collections import Counter
from functools import lru_cache
from math import factorial, prod
import json
from pathlib import Path

MODULUS = 17
A_DATA = [
    ((1,3,13),2), ((1,11,5),1), ((7,8,2),3), ((8,15,11),2),
    ((10,0,7),1), ((10,6,1),1), ((11,10,13),1), ((11,16,7),1),
    ((14,6,14),1), ((14,15,5),1), ((15,1,1),4),
]
B_DATA = [
    ((1,15,1),3), ((7,3,7),2), ((7,16,11),1), ((8,8,1),2),
    ((10,10,14),1), ((10,11,13),1), ((11,1,5),2), ((14,1,2),2),
    ((15,0,2),1), ((15,6,13),2), ((15,8,11),1),
]

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def expand(data):
    return [r for r,m in data for _ in range(m)]

A, B = expand(A_DATA), expand(B_DATA)
require(len(A) == len(B) == 18, "Wrong table sizes")
require(all(sum(r) % MODULUS == 0 for r in A+B), "Nonzero row sum")
require(set(A).isdisjoint(B), "Tables have a common row type")
margin_a = [Counter(r[j] for r in A) for j in range(3)]
margin_b = [Counter(r[j] for r in B) for j in range(3)]
require(margin_a == margin_b, "Incompatible column multisets")
first = tuple(sorted(r[0] for r in A))
bs, cs = tuple(sorted(margin_a[1])), tuple(sorted(margin_a[2]))
c_index = {c:i for i,c in enumerate(cs)}

@lru_cache(maxsize=None)
def count_fillings(position, b_left, c_left):
    if position == len(first):
        return int(not any(b_left) and not any(c_left))
    answer = 0
    a = first[position]
    for i,b in enumerate(bs):
        if not b_left[i]:
            continue
        c = (-a-b) % MODULUS
        j = c_index.get(c)
        if j is None or not c_left[j]:
            continue
        nb, nc = list(b_left), list(c_left)
        nb[i] -= 1
        nc[j] -= 1
        answer += count_fillings(position+1, tuple(nb), tuple(nc))
    return answer

b_start = tuple(margin_a[1][v] for v in bs)
c_start = tuple(margin_a[2][v] for v in cs)
total_fillings = count_fillings(0, b_start, c_start)

def labelings(data):
    answer = 1
    for a, multiplicity in margin_a[0].items():
        denominator = prod(factorial(m) for (x,y,z),m in data if x == a)
        answer *= factorial(multiplicity) // denominator
    return answer

a_fillings, b_fillings = labelings(A_DATA), labelings(B_DATA)
require(a_fillings == 24, "Unexpected A labeling count")
require(b_fillings == 72, "Unexpected B labeling count")
require(total_fillings == a_fillings+b_fillings == 96,
        "An additional compatible table exists")

# A and B account for all 96 fillings. Every compatible table has at least
# one filling when its rows are sorted by their first entries. Hence there
# are exactly two unordered tables. Each has no move of size <= 17.
report = {
    "method": "dynamic_programming_fixed_first_column",
    "labelled_fillings_total": total_fillings,
    "labelled_fillings_A": a_fillings,
    "labelled_fillings_B": b_fillings,
    "dynamic_programming_states": count_fillings.cache_info().currsize,
    "unordered_compatible_tables": 2,
    "no_common_row": True,
    "minimum_nontrivial_move_size": 18,
    "independent_of_rank_calculation": True,
}
print(json.dumps(report, indent=2))
Path(__file__).with_name("audit_Z17_independent_result.json").write_text(
    json.dumps(report, indent=2)+"\n", encoding="utf-8")
