from sympy import *
from Riemann import *

v=symbols("r theta")
u=[v[0]*cos(v[1]), v[0]*sin(v[1])]
g=Matrix([[1, 0],[0, 1]])
g_v=metric_transform(g, u, v)
print(g_v)