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
file_to_write = open('dynamic_answers.txt', 'w')
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
  end_time = time.time()
  from mip import *

  file_to_write.write(filename + '\n')
  file_to_write.write('N = ' + str(N-1) + '\n')
  file_to_write.write('W = ' + str(W) + '\n')
  file_to_write.write('int(W/min_w) = ' + str(int(W / min(w[1:]))) + '\n')
  file_to_write.write('variants = ' + str(len(final_variants)) + '\n')
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