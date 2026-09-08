#!/usr/bin/env python3
"""Optional symbolic audit. Requires SymPy; the other two scripts do not.
Run: python audit_Z17_symbolic.py
"""
from collections import Counter
from itertools import product
import json
from pathlib import Path
import sympy as s
import audit_Z17_independent as data

S = [sorted({r[j] for r in data.A}) for j in range(3)]
E = sorted(r for r in product(*S) if sum(r) % 17 == 0)
M = s.Matrix([[int(r[j] == v) for r in E] for j in range(3) for v in S[j]])
N = M.extract([i for i in range(23) if i not in (15, 22)], list(range(21)))
a,b,d,e,f,g,h = s.symbols('a b d e f g h')
z = s.Matrix([a,b,-a-b,-a,a-d,d,e,-e,f,b-f-g,g,-b,
              g+d,-g,-d,-h-a-b-e,h,a+b+e,-f,
              -g-d+h+a+b+e,-b+f+g-h,-a+d-e])
conditions = [s.expand(t) for t in M*z]
solution = s.linsolve(conditions, [a,b,d,e,f,g,h])
ca, cb = Counter(data.A), Counter(data.B)
kernel = s.Matrix([cb[r]-ca[r] for r in E])
require = data.require
require(len(E) == 22, 'Wrong number of allowed types')
require(M.rank() == 21, 'Wrong rank')
require(M.nullspace() == [kernel], 'Wrong nullspace')
require(N.det() == 17, 'Wrong determinant')
require(solution == s.FiniteSet((2*h,h,-h,-2*h,h,-h,h)),
        'Incorrect hand-written parametrization')
report = {
    'matrix_shape': [23,22],
    'rank_Q': M.rank(),
    'nullspace_basis_B_minus_A': list(map(int, kernel)),
    'minor_det': int(N.det()),
    'symbolic_parametrization': str(solution),
    'engine': 'SymPy exact rational arithmetic',
}
print(json.dumps(report, indent=2))
Path(__file__).with_name('audit_Z17_symbolic_result.json').write_text(
    json.dumps(report, indent=2)+'\n', encoding='utf-8')
