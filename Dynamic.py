import time
def dynamic(A, W, N, w, p, unbounded):
  if unbounded:
    for i in range(0,N):
      for c in range(0,W+1):
        if c<=(w[i]-1):
          A[i][c] = A[i-1][c]
        else:
          A[i][c] = max(A[i-1][c], A[i][c-w[i]]+p[i])
  else:
    for k in range(0,N):
      for s in range(0,W+1):
        if s >= w[k]:
          A[k][s] = max(A[k - 1][s], A[k - 1][s - w[k]] + p[k])
        else:
          A[k][s] = A[k - 1][s]

def findAns(A, k, s, w, unbounded):
  if A[k][s] == 0:
    return answer
  if k>=1:
    if A[k - 1][s] == A[k][s]:
      findAns(A, k - 1, s, w, unbounded)
    else:
      if unbounded:
        findAns(A, k, s - w[k], w, unbounded)
      else:
        findAns(A, k - 1, s - w[k], w, unbounded)
      answer.append(k)
  else:
    return answer

def reverse_answer(N, answer):
  for i in range(0, len(answer)):
    answer[i] = abs(N-answer[i])
  return answer

import glob
filenames1 = glob.glob("./tests/test?.txt")
filenames2 = glob.glob("./tests/test??.txt")
filenames = filenames1 + filenames2

for filename in filenames:
  start_time = time.time()
  f = open(filename)
  N = int(f.readline())
  W = int(f.readline())
  w = [0]
  b = []
  for i in range(N):
    new_line = f.readline().split(' ')
    w.append(int(new_line[0]))
    b.append(int(new_line[1]))

  N = N + 1
  p = [1 for i in range(N)]
  p[0] = 0
  A = [[0 for j in range(W+1)] for i in range(0,N)]

  min_w = min(w[1:N])

  dynamic(A, W, N, w, p, unbounded = False)
  answers = []
  for i in range(1,N):
    for j in range(min_w ,W+1):
      answer = []
      findAns(A, i, j, w, unbounded=False)
      answer.sort()
      if (len(answer) > 0) and (answer not in answers):
        answers.append(answer)

  A = [[0 for j in range(W + 1)] for i in range(0, N)]
  dynamic(A, W, N, w, p, unbounded = True)
  for i in range(1,N):
    for j in range(min_w ,W+1):
      answer = []
      findAns(A, i, j, w, unbounded=True)
      answer.sort()
      if (len(answer) > 0) and (answer not in answers):
        answers.append(answer)

  w = w[1:N]
  w.reverse()
  w = [0] + w

  A = [[0 for j in range(W + 1)] for i in range(0, N)]
  dynamic(A, W, N, w, p, unbounded = False)
  for i in range(1, N):
    for j in range(min_w, W + 1):
      answer = []
      findAns(A, i, j, w, unbounded=False)
      answer = reverse_answer(N, answer)
      answer.sort()
      if (len(answer) > 0) and (answer not in answers):
        answers.append(answer)

  A = [[0 for j in range(W + 1)] for i in range(0, N)]
  dynamic(A, W, N, w, p, unbounded = True)
  for i in range(1, N):
    for j in range(min_w, W + 1):
      answer = []
      findAns(A, i, j, w, unbounded=True)
      answer.sort()
      answer = reverse_answer(N, answer)
      if (len(answer) > 0) and (answer not in answers):
        answers.append(answer)


  variants = []
  for answer in answers:
    variant = [0 for j in range(N-1)]
    for j in range(0, len(answer)):
      variant[answer[j]-1] += 1
    variants.append(variant)


  sum_of_columns = [sum(row[i] for row in variants) for i in range(len(variants[0]))]
  for i in range(len(sum_of_columns)):
    if sum_of_columns[i] == 0:
      new_variant = [0 for j in range(len(sum_of_columns))]
      new_variant[i] = int(W/w[i])
      variants.append(new_variant)

  c2 = []
  final_variants = []
  for i in variants:
    weight = 0
    for j in range(0,len(i)):
      weight += i[j]*w[j+1]
    if (W-weight)>=0:
      final_variants.append(i)
      c2.append(W-weight)

  time_for_searching_variants = time.time() - start_time
  print(filename)
  N = N - 1
  print('N = ', N)
  print('W = ', W)
  print('int(W/min_w) = ', int(W/min_w))
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