import time
import random
random.seed(0)
def fit_fun(population, w, num_items, W):
    fit_population = []
    value_population = []
    for chromosomes in population:
        total_value = sum([chromosomes[i] * w[i] for i in range(0, num_items)])
        while (total_value > W):
            indexes = [i for i, v in enumerate(chromosomes) if v > 0]
            ind_for_change = random.choices(indexes)[0]
            chromosomes[ind_for_change] -= 1
            total_value = sum([chromosomes[i] * w[i] for i in range(num_items)])
        value_population.append(total_value)
        fit_population.append(sum([chromosomes[i] * w[i] for i in range(num_items)]))
    return (value_population, fit_population)

def check_population(value_population, fit_population, num_items):
    max_value = max(fit_population)
    max_ind = fit_population.index(max_value)
    max_count = fit_population.count(max_value)
    if (max_count >= (0.6 * num_items)):
        return True
    else:
        return False

def single_point_crossover(chromosome_x, chromosome_y):
    crossover_ind = random.randint(0, len(chromosome_x) - 1)
    return chromosome_x[:crossover_ind] + chromosome_y[crossover_ind:]

def group_selection(fit_population, generation_size):
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

def mutation(new_population, w):
    for chromosomes in new_population:
        choice = random.choices([0, 1], weights=[0.9, 0.1], k=1)
        if (choice == 1):
            random_chromosome = random.randint(0, int(W/min(w)))
            ind = random.randint(0, len(chromosomes) - 1)
            while chromosomes[ind] == random_chromosome:
                random_chromosome = random.randint(0, int(W / min(w)))
            chromosomes[ind] = random_chromosome
    return new_population


import glob
filenames1 = glob.glob("./tests/test?.txt")
filenames2 = glob.glob("./tests/test??.txt")
filenames = filenames1 + filenames2
file_to_write = open('genetic_answers.txt', 'w')
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
    c2 = []
    p = [1 for i in range(0, len(w))]
    num_chromosomes = N ** 3
    start_num_chromosomes = N ** 3
    num_items = len(w)
    variants = []
    population = [random.choices([i for i in range(0, int(W/min(w))+1)], k=num_items) for i in range(0, start_num_chromosomes)]
    while (1):
        value_population, fit_population = fit_fun(population, w, num_items, W)
        flag = check_population(value_population, fit_population, num_items)
        if flag:
            variants_from_genetic = population
            break
        else:
            new_population = group_selection(fit_population, num_chromosomes)
            population = mutation(new_population, w)

    variants = []
    for variant in variants_from_genetic:
        if variant not in variants:
            variants.append(variant)

    sum_of_columns = [sum(row[i] for row in variants) for i in range(len(variants[0]))]
    for i in range(len(sum_of_columns)):
        if sum_of_columns[i] == 0:
            new_variant = [0 for j in range(len(sum_of_columns))]
            new_variant[i] = int(W / w[i])
            variants.append(new_variant)

    final_variants = []
    optimal_variants = 0
    for variant in variants:
        if variant not in final_variants:
            final_variants.append(variant)
            weight = 0
            for j in range(0, len(variant)):
                weight += w[j] * variant[j]
            c2.append(W - weight)
            if (W-weight)<min(w):
                optimal_variants+=1

    end_time = time.time()
    from mip import *

    file_to_write.write(filename + '\n')
    file_to_write.write('N = ' + str(N - 1) + '\n')
    file_to_write.write('W = ' + str(W) + '\n')
    file_to_write.write('int(W/min_w) = ' + str(int(W / min(w))) + '\n')
    file_to_write.write('variants = ' + str(len(final_variants)) + '\n')
    file_to_write.write('optimal variants = ' + str(optimal_variants) + '\n')
    file_to_write.write('time for searching variants = ' + str(end_time - start_time) + '\n')
    m = Model()
    a = np.array(final_variants).transpose()
    c1 = [1 for i in range(a.shape[1])]
    I = range(0, len(final_variants))
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
    file_to_write.write('min_of_waste = ' + str(m.objective_value) + '\n\n\n')
file_to_write.close()
