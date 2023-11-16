from sympy import *
from Riemann import *


u_space=symbols("r theta phi")
# r 0 || theta 1 || phi 2


u=symbols("r theta phi t")

k=symbols("k")
g_space=Matrix([[1/(1-k/u[0]), 0, 0],[0, u[0]**2 ,0], [0,0,u[0]**2*(sin(u[1]))**2]])
g=Matrix([[1/(1-k/u[0]), 0, 0, 0],
        [0, u[0]**2 ,0, 0],
        [0,0,u[0]**2*(sin(u[1]))**2, 0],
        [0,0,0,-(1-k/u[0])]])
R=Riemann_Tensor(g_space, u_space)
#Rc=Ricci_Tensor(g, u)
#R=Ricci_Scalar(g, u)
print(R)