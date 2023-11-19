from sympy import *
from Riemann import *

v=symbols("p q")
u=[v[0]*v[1], v[0]+v[1]] #x, y
g=Matrix([[1, 0],[0, 1]])
g_v=metric_transform(g, u, v)
print(g_v)

R=Riemann_Tensor(g_v, v)
print(R)