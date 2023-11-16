from sympy import *
from Riemann import *


N=2
u=symbols("x y")
r=symbols("r")
g=Matrix([[(r**2)*(sin(u[1]))**2, 0],[0, r**2]])
#Chris=Christoffel_Tensor(g, u)
R=Ricci_Scalar(g, u)
print(R)