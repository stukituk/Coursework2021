import time
import glob

filenames1 = glob.glob("./tests/test?.txt")
filenames2 = glob.glob("./tests/test??.txt")
filenames = filenames1 + filenames2
file_to_write = open('brute_force_answers.txt', 'w')
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
    import numpy as np

    end_time = time.time()
    from mip import *
    file_to_write.write(filename + '\n')
    file_to_write.write('N = ' + str(N) + '\n')
    file_to_write.write('W = ' + str(W) + '\n')
    file_to_write.write('int(W/min_w) = ' +  str(int(W/min(w))) + '\n')
    file_to_write.write('variants = ' + str(len(final_variants)) + '\n')
    file_to_write.write('time for searching variants = ' + str(end_time-start_time) + '\n')
    m = Model()
    a = np.array(final_variants).transpose()
    c1 = [1 for i in range(a.shape[1])]
    I = range(0,len(final_variants))
    b = np.array(b)
    b_new = 0
    x = [m.add_var(var_type=INTEGER) for i in I]

    for A in a:
        m += xsum(A[i] * x[i] for i in I) >= b[b_new]
        b_new += 1
    m.objective = minimize(xsum(c1[i] * x[i] for i in I))
    m.optimize(max_seconds=60)
    file_to_write.write('min_of_sheets = ' + str(m.objective_value) + '\n')
    m.objective = minimize(xsum(c2[i] * x[i] for i in I))
    m.optimize(max_seconds=60)
    file_to_write.write('min_of_waste = ' + str(m.objective_value) + '\n\n')
file_to_write.close()

