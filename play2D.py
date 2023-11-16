from sympy import *
from Riemann import *


N=2
u=symbols("r theta")
g=Matrix([[1, 0],[0, u[0]**4]])
#Chris=Christoffel_Tensor(g, u)

R=Riemann_Tensor(g, u)
print(R)