#!/usr/bin/env python3
"""Точная проверка контрпримера для G = Z/17Z.

Запуск: python verify_phylogenetic_Z17.py
Требуется только стандартная библиотека Python 3.9+.

Проверяются независимо:
  1) линейные уравнения для кратностей строк над рациональными числами;
  2) полный перебор всех таблиц с заданными наборами элементов в столбцах.

Равенства кратностей — обычные целочисленные равенства, НЕ равенства
по модулю 17. По модулю 17 проверяются только суммы элементов строк.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path
from typing import Sequence

Q = 17
A = [
    (1, 3, 13), (1, 3, 13), (1, 11, 5),
    (7, 8, 2), (7, 8, 2), (7, 8, 2),
    (8, 15, 11), (8, 15, 11),
    (10, 0, 7), (10, 6, 1),
    (11, 10, 13), (11, 16, 7),
    (14, 6, 14), (14, 15, 5),
    (15, 1, 1), (15, 1, 1), (15, 1, 1), (15, 1, 1),
]
B = [
    (1, 15, 1), (1, 15, 1), (1, 15, 1),
    (7, 3, 7), (7, 3, 7), (7, 16, 11),
    (8, 8, 1), (8, 8, 1),
    (10, 10, 14), (10, 11, 13),
    (11, 1, 5), (11, 1, 5),
    (14, 1, 2), (14, 1, 2),
    (15, 0, 2), (15, 6, 13), (15, 6, 13), (15, 8, 11),
]


def require(condition: bool, message: str) -> None:
    """В отличие от assert, проверка не отключается флагом Python -O."""
    if not condition:
        raise RuntimeError(message)


def rref(matrix: Sequence[Sequence[int]]) -> tuple[list[list[Fraction]], list[int]]:
    """Приведенная ступенчатая форма, исключительно точные дроби."""
    result = [[Fraction(value) for value in row] for row in matrix]
    row_count, col_count = len(result), len(result[0])
    pivots: list[int] = []
    lead_row = 0
    for col in range(col_count):
        pivot = next((i for i in range(lead_row, row_count) if result[i][col]), None)
        if pivot is None:
            continue
        result[lead_row], result[pivot] = result[pivot], result[lead_row]
        value = result[lead_row][col]
        result[lead_row] = [entry / value for entry in result[lead_row]]
        for i in range(row_count):
            if i != lead_row and result[i][col]:
                factor = result[i][col]
                result[i] = [
                    left - factor * right
                    for left, right in zip(result[i], result[lead_row])
                ]
        pivots.append(col)
        lead_row += 1
        if lead_row == row_count:
            break
    return result, pivots


def determinant(matrix: Sequence[Sequence[int]]) -> Fraction:
    """Определитель квадратной матрицы с точной рациональной арифметикой."""
    result = [[Fraction(value) for value in row] for row in matrix]
    n = len(result)
    require(all(len(row) == n for row in result), 'Матрица не квадратная.')
    answer = Fraction(1)
    for col in range(n):
        pivot = next((i for i in range(col, n) if result[i][col]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            result[col], result[pivot] = result[pivot], result[col]
            answer *= -1
        diagonal = result[col][col]
        answer *= diagonal
        for i in range(col + 1, n):
            factor = result[i][col] / diagonal
            for j in range(col + 1, n):
                result[i][j] -= factor * result[col][j]
            result[i][col] = Fraction(0)
    return answer


def enumerate_tables(
    allowed: list[tuple[int, int, int]],
    margins: list[Counter[int]],
) -> tuple[set[tuple[int, ...]], int]:
    """Полный перебор кратностей; не использует ранг или вектор ядра.

    Для каждого значения первого столбца перебираются все распределения
    его кратности между разрешенными строками. Ветвь отбрасывается только
    при превышении кратности в столбце 2 или 3. Таким образом, ни одна
    совместимая таблица не пропускается.
    """
    first_values = sorted(margins[0])
    groups = [[i for i, row in enumerate(allowed) if row[0] == a]
              for a in first_values]
    remaining_second = margins[1].copy()
    remaining_third = margins[2].copy()
    counts = [0] * len(allowed)
    solutions: set[tuple[int, ...]] = set()
    visited = 0

    def assign_group(group_index: int) -> None:
        nonlocal visited
        visited += 1
        if group_index == len(groups):
            if (all(v == 0 for v in remaining_second.values()) and
                    all(v == 0 for v in remaining_third.values())):
                solutions.add(tuple(counts))
            return
        indices = groups[group_index]
        total = margins[0][first_values[group_index]]

        def assign_row(position: int, remaining_total: int) -> None:
            if position == len(indices):
                if remaining_total == 0:
                    assign_group(group_index + 1)
                return
            index = indices[position]
            _, b, c = allowed[index]
            low = remaining_total if position == len(indices) - 1 else 0
            high = min(remaining_total, remaining_second[b], remaining_third[c])
            for value in range(low, high + 1):
                counts[index] = value
                remaining_second[b] -= value
                remaining_third[c] -= value
                assign_row(position + 1, remaining_total - value)
                remaining_second[b] += value
                remaining_third[c] += value
            counts[index] = 0

        assign_row(0, total)

    assign_group(0)
    return solutions, visited


def verify(method: str = 'all') -> dict:
    require(len(A) == len(B) == 18, 'Ожидались две таблицы по 18 строк.')
    require(all(len(row) == 3 for row in A + B), 'Неверная ширина таблицы.')
    require(all(0 <= g < Q for row in A + B for g in row), 'Неверный остаток.')
    require(all(sum(row) % Q == 0 for row in A + B), 'Есть ненулевая сумма строки.')
    margins = [Counter(row[j] for row in A) for j in range(3)]
    require(margins == [Counter(row[j] for row in B) for j in range(3)],
            'Наборы элементов соответствующих столбцов различаются.')
    require(not (set(A) & set(B)), 'Обнаружена общая строка.')

    # Проверяются все 17^3 тройки, а не только заранее выписанные 22.
    allowed = [row for row in product(range(Q), repeat=3)
               if sum(row) % Q == 0
               and all(margins[j][row[j]] > 0 for j in range(3))]
    require(len(allowed) == 22, 'Число разрешенных типов строк не равно 22.')
    counts_a = tuple(Counter(A)[row] for row in allowed)
    counts_b = tuple(Counter(B)[row] for row in allowed)
    kernel = tuple(a - b for a, b in zip(counts_a, counts_b))
    result = {'group_order': Q, 'rows_per_table': len(A),
              'columns': 3, 'allowed_row_types': len(allowed),
              'compatible_column_multisets': True,
              'zero_row_sums': True, 'no_common_row_type': True,
              'kernel_A_minus_B': kernel}

    print('Группа: Z/17Z. Размер каждой таблицы: 18 x 3.')
    print('Суммы строк нулевые; наборы элементов столбцов совпадают.')
    print('Общих типов строк у A и B нет.')
    print('Полная проверка 17^3 троек: разрешено ровно 22 типа строк.')
    for j, margin in enumerate(margins, 1):
        print(f'Столбец {j}: {sorted(margin.items())}')

    if method in ('all', 'linear'):
        labels = [(j, value) for j in range(3) for value in sorted(margins[j])]
        matrix = [[int(row[j] == value) for row in allowed] for j, value in labels]
        reduced, pivots = rref(matrix)
        require(len(pivots) == 21, 'Ранг не равен 21.')
        require(pivots == list(range(21)), 'Неожиданные ведущие столбцы.')
        require(all(sum(x * y for x, y in zip(row, kernel)) == 0 for row in matrix),
                'Вектор A-B не принадлежит ядру.')
        require(kernel[-1] == -1, 'Последняя координата ядра не равна -1.')
        require(all(reduced[i][-1] == kernel[i] for i in range(21)),
                'Приведенная матрица не совпадает с сертификатом ядра.')
        # Удаляются строки 16 и 23 (нумерация с единицы), столбец 22.
        selected_rows = [i for i in range(23) if i not in (15, 22)]
        minor = [[matrix[i][j] for j in range(21)] for i in selected_rows]
        det = determinant(minor)
        require(det == 17, 'Определитель указанного минора не равен 17.')
        # x = A - t*(A-B). Поскольку x_22=t, параметр целый.
        # x_2=1-t и x_22=t неотрицательны, поэтому t=0 или t=1.
        require(counts_a[-1] == 0 and counts_b[-1] == 1, 'Неожиданная строка 22.')
        require(counts_a[1] == 1 and counts_b[1] == 0, 'Неожиданная строка 2.')
        for t, expected in [(0, counts_a), (1, counts_b)]:
            require(tuple(a - t * u for a, u in zip(counts_a, kernel)) == expected,
                    'Неверная параметризация решений.')
        result.update(exact_rank=21, nullity=1, certified_minor_determinant=17,
                      fiber_size_by_linear_algebra=2)
        print('Точная рациональная арифметика: ранг 21, размерность ядра 1.')
        print('Указанный минор размера 21 x 21 имеет определитель 17.')
        print('Линейный вывод: x=(1-t)A+tB, целое t; 0<=t<=1. Ровно два решения.')

    if method in ('all', 'enumerate'):
        solutions, visited = enumerate_tables(allowed, margins)
        require(solutions == {counts_a, counts_b},
                'Полный перебор не подтвердил, что решения — только A и B.')
        result.update(fiber_size_by_enumeration=len(solutions),
                      enumeration_group_states=visited)
        print(f'Независимый полный перебор: ровно {len(solutions)} таблицы — A и B.')
        print(f'Посещено состояний между группами: {visited}.')

    print('Вывод: нет перехода из A при замене не более 17 из 18 строк.')
    print('Следовательно, phi(Z/17Z) >= 18 > 17.')
    print('Все проверки пройдены.')
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--method', choices=('all', 'linear', 'enumerate'), default='all',
                        help='Метод проверки; по умолчанию оба независимо.')
    parser.add_argument('--report', type=Path, help='Путь для необязательного отчета JSON.')
    args = parser.parse_args()
    result = verify(args.method)
    if args.report is not None:
        args.report.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n',
                               encoding='utf-8')


if __name__ == '__main__':
    main()
