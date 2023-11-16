from sympy import *
init_printing(use_unicode=True)


def Christoffel(i, j, k, g, u):
    N=len(u)
    g_inv=g.inv()
    result=0
    for t1 in range(N):
        result+=g_inv[k, t1]*(diff(g[t1, i], u[j])+diff(g[j, t1], u[i])-diff(g[i,j], u[t1]))
    return expand(result/2)

def Christoffel_Tensor(g, u):
    N=len(u)
    return [[[Christoffel(i, j, k, g, u) for k in range(N)] for j in range(N)] for i in range(N)]

def Riemann(c, a, b, d, g, u):
    N=len(u)
    Chris=Christoffel_Tensor(g, u)
    result=diff(Chris[b][c][d], u[a])-diff(Chris[a][c][d], u[b])
    for t1 in range(N):
        result+=Chris[b][c][t1]*Chris[a][t1][d]-Chris[a][c][t1]*Chris[b][t1][d]
    return result

def Riemann_2(c, a, b, d, g, u, Chris):
    N=len(u)
    result=diff(Chris[b][c][d], u[a])-diff(Chris[a][c][d], u[b])
    for t1 in range(N):
        result+=Chris[b][c][t1]*Chris[a][t1][d]-Chris[a][c][t1]*Chris[b][t1][d]
    return result

def Riemann_Tensor(g, u):
    N=len(u)
    Chris=Christoffel_Tensor(g, u)
    return [[[[Riemann_2(l,i, j, k, g, u, Chris) for k in range(N)] for j in range(N)] for i in range(N)] for l in range(N)]

def Ricci(a, b, g, u):
    N=len(u)
    R=Riemann_Tensor(g, u)
    result=0
    for rho in range(N):
        result+=R[a][rho][b][rho]
    return result 

def Ricci_2(a, b, g, u, R):
    N=len(u)
    result=0
    for rho in range(N):
        result+=R[a][rho][b][rho]
    return result

def Ricci_Tensor(g, u):
    N=len(u)
    R=Riemann_Tensor(g, u)
    return [[Ricci_2(a, b, g, u, R) for b in range(N)] for a in range(N)] 

def Ricci_Scalar(g, u):
    g_inv=g.inv()
    N=len(u)
    result=0 
    R=Ricci_Tensor(g, u)
    for a in range(N):
        for b in range(N):
            result+=R[a][b]*g_inv[a,b]
    return result

N=2
u=[0,0]
u=symbols("x y")
r=symbols("r")
g=Matrix([[(r**2)*(sin(u[1]))**2, 0],[0, r**2]])
#Chris=Christoffel_Tensor(g, u)
R=Ricci_Scalar(g, u)
print(R)