from cmath import pi
from math import tan, sqrt
import numpy as np
# w = 1 /125
# w = 2 * 2 * tan(pi * w / 2)

# a = 16 + 4*2 ** 0.5 * w + w ** 2
# b = -32 + 2 * w ** 2
# c = 16 - 4*2 ** 0.5 * w + w ** 2

# print(np.roots([a, b, c]))
# print(c/a)


def calc(w1, w2):
    f = 250
    w = tan((w2 - w1) * pi / f)
    omega_c = tan(sqrt(w1 * w2) * pi / f)

    b = 1 + 2 ** 0.5 * w + 2 * omega_c ** 2 + w ** 2 + 2 ** 0.5 * omega_c ** 2 * w + omega_c ** 4
    a = w ** 2
    a = 6 - 2 * (2 * omega_c * omega_c + w * w) + 6 * omega_c ** 4
    a = 1 - 2 ** 0.5 * w + 2 * omega_c ** 2 + w ** 2 - 2 ** 0.5 * omega_c ** 2 * w + omega_c ** 4

    print(a / b)
    print(b, a)

w1 = 8
w2 = 20

calc(w1, w2)