import time
import random
random.seed(0)
def fit_fun(population, w, num_items, W):
    fit_population = []
    value_population = []
    for chromosomes in population:
        total_value = sum([chromosomes[i] * w[i] for i in range(0, num_items)])
        while (total_value > W):
            indexes = [i for i, v in enumerate(chromosomes) if v == 1]
            ind_for_change = random.choices(indexes)[0]
            chromosomes[ind_for_change] = 0
            total_value = sum([chromosomes[i] * w[i] for i in range(num_items)])
        value_population.append(total_value)
        fit_population.append(sum([chromosomes[i] * w[i] for i in range(num_items)]))
    return (value_population, fit_population)

def check_population(value_population, fit_population, num_items):
    max_value = max(fit_population)
    max_ind = fit_population.index(max_value)
    max_count = fit_population.count(max_value)
    if (max_count >= (0.6 * num_items)):
        return (True, max_ind)
    else:
        return (False, max_ind)

def single_point_crossover(chromosome_x, chromosome_y):
    crossover_ind = random.randint(0, len(chromosome_x) - 1)
    return chromosome_x[:crossover_ind] + chromosome_y[crossover_ind:]

def group_selection(it_population, generation_size):
    sorted_population = sorted(range(len(fit_population)), key=lambda k: fit_population[k], reverse=True)
    group_size = int(len(sorted_population) // 4)
    random_weights = [0.5] * group_size + [0.3] * group_size + [0.15] * group_size + [0.05] * group_size
    diff = len(sorted_population) - len(random_weights)
    for i in range(diff):
        random_weights.append(0.05)
    new_population = []
    for i in range(generation_size):
        x, y = random.choices(sorted_population, weights=random_weights, k=2)
        new_chromosome = single_point_crossover(population[x], population[y])
        new_population.append(new_chromosome)
    return new_population

def mutation(new_population):
    for chromosomes in new_population:
        choice = random.choices([0, 1], weights=[0.9, 0.1], k=1)
        if (choice == 1):
            ind = random.randint(0, len(chromosomes) - 1)
            if (chromosomes[ind] == 0):
                chromosomes[ind] = 1
            else:
                chromosomes[ind] = 0
    return new_population

import glob
filenames1 = glob.glob("./tests/test?.txt")
filenames2 = glob.glob("./tests/test??.txt")
filenames = filenames1 + filenames2

for filename in filenames:
    start_time = time.time()
    f = open(filename)
    N = int(f.readline())
    W = int(f.readline())
    b = []
    w = []
    for i in range(N):
        new_line = f.readline().split(' ')
        w.append(int(new_line[0]))
        b.append(int(new_line[1]))
    w_for_unbounded = []
    indexses = []
    c2 = []
    for i in range(0, len(w)):
        count = int(W / w[i])
        for j in range(0, count):
            indexses.append(i)
            w_for_unbounded.append(w[i])
    p = [1 for i in range(0, len(w_for_unbounded))]
    start_num_chromosomes = len(p)
    if len(p)>=100:
        num_chromosomes = len(p) * 10
    else:
        num_chromosomes = len(p) ** 2
    num_items = len(w_for_unbounded)
    variants = []
    population = [random.choices([0, 1], k=num_items) for i in range(0, num_chromosomes)]
    while (1):
        value_population, fit_population = fit_fun(population, w_for_unbounded, num_items, W)
        flag, max_ind = check_population(value_population, fit_population, num_items)
        if flag:
            variants_from_genetic = population
            break
        else:
            new_population = group_selection(fit_population, num_chromosomes)
            population = mutation(new_population)

    variants = []
    for variant in variants_from_genetic:
        new_var = [0 for i in range(0, len(w))]
        for j in range(0, len(variant)):
            if variant[j] == 1:
                new_var[indexses[j]] += 1
        if new_var not in variants:
            variants.append(new_var)

    sum_of_columns = [sum(row[i] for row in variants) for i in range(len(variants[0]))]
    for i in range(len(sum_of_columns)):
        if sum_of_columns[i] == 0:
            new_variant = [0 for j in range(len(sum_of_columns))]
            new_variant[i] = int(W / w[i])
            variants.append(new_variant)

    final_variants = []
    for variant in variants:
        if variant not in final_variants:
            final_variants.append(variant)
            weight = 0
            for j in range(0, len(variant)):
                weight += w[j] * variant[j]
            c2.append(W - weight)

    time_for_searching_variants = time.time() - start_time
    print(filename)
    N = N
    print('N = ', N)
    print('W = ', W)
    print('int(W/min_w) = ', int(W/min(w)))
    print('time for searching variants = ', time_for_searching_variants)
    print('variants = ', len(final_variants))
    from scipy.optimize import linprog
    import numpy as np

    a = - np.array(final_variants).transpose()

    c1 = [1 for i in range(a.shape[1])]

    b = - np.array(b)

    bnd = [(0, float("inf"))]
    start_time = time.time()
    res1 = linprog(c1, A_ub=a, b_ub=b, bounds=bnd, method='simplex')
    import math

    x1 = res1.x
    for i in range(0, len(x1)):
        x1[i] = math.ceil(x1[i])
    time_for_searching_the_min_of_sheets = time.time() - start_time
    print('time for searching the min of sheets = ', time_for_searching_the_min_of_sheets)
    print('min of sheets = ', sum(x1))
    import math

    start_time = time.time()
    res2 = linprog(c2, A_ub=a, b_ub=b, bounds=bnd, method='simplex')
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