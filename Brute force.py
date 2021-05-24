import time
import glob

filenames1 = glob.glob("./tests/test?.txt")
filenames2 = glob.glob("./tests/test??.txt")
filenames = filenames1 + filenames2
for filename in filenames:
    start_time = time.time()
    f = open(filename)
    N = int(f.readline())
    W = int(f.readline())
    w = []
    b = []
    c2 = []
    optimal_variants = []
    for i in range(N):
        new_line = f.readline().split(' ')
        w.append(int(new_line[0]))
        b.append(int(new_line[1]))
    indexses = range(0, len(w))
    import itertools

    variants = []
    for i in range(1, int(W / min(w)) + 1):
        variants = itertools.combinations_with_replacement(indexses, i)
        for q in variants:
            weight = 0
            for j in range(0, len(q)):
                weight += w[q[j]]
            if (weight <= W) and (weight + min(w) > W):
                optimal_variants.append(q)
                c2.append(W - weight)

    final_variants = []
    for i in optimal_variants:
        variant = [0 for j in range(0, len(w))]
        for j in range(0, len(i)):
            variant[i[j]] += 1
        final_variants.append(variant)
    time_for_searching_variants = time.time() - start_time
    print(filename)
    print('N = ', N)
    print('W = ', W)
    print('int(W/min_w) = ', int(W/min(w)))
    print('time for searching variants = ', time_for_searching_variants)
    print('variants = ', len(final_variants))
    from scipy.optimize import linprog
    import numpy as np

    A = - np.array(final_variants).transpose()

    c1 = [1 for i in range(A.shape[1])]

    b = - np.array(b)

    bnd = [(0, float("inf"))]
    start_time = time.time()
    res1 = linprog(c1, A_ub=A, b_ub=b, bounds=bnd, method='simplex')
    import math

    x1 = res1.x
    for i in range(0, len(x1)):
        x1[i] = math.ceil(x1[i])
    time_for_searching_the_min_of_sheets = time.time() - start_time
    print('time for searching the min of sheets = ', time_for_searching_the_min_of_sheets)
    print('min of sheets = ', sum(x1))
    import math
    start_time = time.time()
    res2 = linprog(c2, A_ub=A, b_ub=b, bounds=bnd, method='simplex')
    x2 = res2.x
    for i in range(0, len(x2)):
        x2[i] = math.ceil(x2[i]) * c2[i]
    time_for_searching_the_min_of_waste = time.time() - start_time
    print('time for searching the min of waste = ', time_for_searching_the_min_of_waste)
    print('min of waste = ', sum(x2))
    end_time = time.time() - start_time
    print('total_time = ',
          time_for_searching_variants + time_for_searching_the_min_of_waste + time_for_searching_the_min_of_sheets)
    print()


