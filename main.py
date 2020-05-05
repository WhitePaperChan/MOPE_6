import random
import numpy
import copy
from scipy.stats import t
from scipy.stats import f
import math

def cochran(f1, f2, q):
    fish = f.isf(q/f2, f1, (f2 - 1)*f1)
    result = fish/(fish + f2 - 1)
    return result

def print_line(m):
    print("-" * 12 * (m+1))

def coeffs_criterias(yi, x, y):
    global b, m
    k = len(x[0])
    mx = [[] for i in range(len(x) + 1)]
    mx[0].append(k)
    for i in range(1, len(x) + 1):
        suma = round(sum(x[i-1]), 5)
        mx[0].append(suma)
        mx[i].append(suma)
        for j in range(0, len(x)):
            mx[i].append(round(sum([round(x[i-1][l] * x[j][l], 5) for l in range(k)]), 5))

    det = numpy.linalg.det(mx)
    delta = round(det, 5)

    my = [round(sum(yi), 5)]
    for i in range(len(x)):
         my.append(round(sum([yi[j]*x[i][j] for j in range(k)]), 5))

    b = [copy.deepcopy(mx) for i in range(len(x) + 1)]
    for i in range(len(x) + 1):
         for j in range(len(x) + 1):
             b[i][j][i] = my[j]
         b[i] = round(numpy.linalg.det(b[i])/delta,5)
         print("b" + str(i) + ": " + str(b[i]))



    S2 = []
    for i in range(len(y)):
        S2.append(sum([(y[i][j] - yi[i])**2 for j in range(len(y[i]))]))
        S2[i] = round(S2[i]/len(y[i]), 3)
    print("S2: " + str(S2))

    Gp = round(max(S2)/sum(S2), 3)
    print("Gp: " + str(Gp))

    f1 = m - 1
    f2 = k

    print("f1:" + str(f1))
    print("f2:" + str(f2))

    alpha = 0.05

    Gcr = round(cochran(f1, f2, alpha), 4)
    print("Gcr: " + str(Gcr))
    if Gp < Gcr:
        print("Cochran's C: OK")
    else:
        print("Cochran's C: :(")
        m += 1
        return generate_y(x)

    S2v = sum(S2) / 4

    S2b = round(S2v / (4 * m), 3)
    Sb = round(math.sqrt(S2b), 3)

    f3 = f1 * f2
    print("f3: " + str(f3))
    tcr = round(t.ppf(1 - alpha / 2, df=f3), 3)
    print("t: " + str(tcr))
    bs = []
    ts = []
    d = 0
    bs.append(round(sum([yi[j] for j in range(len(yi))]) / len(yi), 3))

    ts.append(round(bs[0] / Sb, 3))
    if ts[0] < 0:
        ts[0] *= -1
    if ts[0] > tcr:
        ts[0] = True
        d += 1
    else:
        ts[0] = False
    for i in range(len(x)):
        bs.append(round(sum([yi[j] * x[i][j] for j in range(len(yi))]) / len(yi), 3))
        ts.append(round(bs[i+1] / Sb, 3))
        if ts[i+1] < 0:
            ts[i+1] *= -1
        if ts[i+1] > tcr:
            ts[i+1] = True
            d += 1
        else:
            ts[i+1] = False

    print("Чи значимі b: " + str(ts))

    f4 = k - d
    print("f4: " + str(f4))
    yj = []
    b0 = []
    for i in range(len(b)):
        if ts[i]:
            b0.append(b[i])
        else:
            b0.append(0)
    for j in range(k):
        yj.append(round(b0[0] + sum([x[i-1][j] * b0[i] for i in range(1, len(b0))]), 3))
    print("yj: " + str(yj))

    S2ad = round(m * sum([(yj[i] - yi[i]) ** 2 for i in range(4)]) / f4, 3)
    Fp = round(S2ad / S2v, 3)
    print("Fp: " + str(Fp))
    Fcr = round(f.ppf(1 - alpha, f4, f3), 1)
    print("Fcr: " + str(Fcr))
    if Fp < Fcr:
        print("F-criteria: OK")
    else:
        print("F-criteria: :(")
        start(x)

def generate_y(x):
    print_line(m)
    print("| " + '{:<10}'.format(""), end="")
    for i in range(1, m + 1):
        print("| " + '{:<10}'.format("yi" + str(i)), end="")
    print("|")
    print_line(m)

    y = []
    k = len(x[0])
    for j in range(1, k + 1):
        y.append([])
        print("| " + '{:<10}'.format(j), end="")
        for i in range(1, m + 1):
            r = round(8 + 5.3 * x[0][j - 1] +
                      0.5 * x[1][j - 1] +
                      5.6 * x[2][j - 1] +
                      3.2 * x[7][j - 1] +
                      0.7 * x[8][j - 1] +
                      4.1 * x[9][j - 1] +
                      8.9 * x[3][j - 1] +
                      0.5 * x[4][j - 1] +
                      1.5 * x[5][j - 1] +
                      1.2 * x[6][j - 1] +
                      random.random() * 10 - 5, 3)
            y[j - 1].append(r)
            print("| " + '{:<10}'.format(r), end="")
        print("|")
        print_line(m)

    yi = []
    for i in range(k):
        yi.append(round(1 / m * sum(y[i]), 3))

    print("y (середні): " + str(yi))

    coeffs_criterias(yi, x, y)

def start(x):
    global m
    while True:
        m = input("m (integer) or cancel:")
        if m.isnumeric():
            print("OK")
            m = int(m)
            break
        elif m == "cancel":
            return
        else:
            print("m must be integer")

    generate_y(x)


x = [
    [-30, -30, -30, -30, 0, 0, 0, 0, -40.95, 10.95, -15, -15, -15, -15],
    [10, 10, 60, 60, 10, 10, 60, 60, 35, 35, -8.25, 78.25, 35, 35],
    [10, 35, 10, 35, 10, 35, 10, 35, 22.5, 22.5, 22.5, 22.5, 0.875, 44.125]
]

k = len(x[0])
x.append([round(x[0][i] * x[1][i], 6) for i in range(k)])
x.append([round(x[0][i] * x[2][i], 6) for i in range(k)])
x.append([round(x[1][i] * x[2][i], 6) for i in range(k)])
x.append([round(x[0][i] * x[1][i] * x[2][i], 6) for i in range(k)])
for j in range(3):
    x.append([round(x[j][i] ** 2, 6) for i in range(k)])
start(x)