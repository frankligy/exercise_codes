from gekko import GEKKO
import numpy as np
import math

'''
Let's set:
DNAdamage as x,
Mdm2 cyto as y,
Mdm2 nuclear as z,
p53 as p
'''

# all the parameters
k_prime_d53 = 0.27
k_prime2_d53 = 8.25
theta = 0.83
k_prime_d2 = 0.05
k_prime_s53 = 0.6
k_prime2_s53 = 2.56
J_s53 = 0.45
k_prime_s2 = 0.15
k_prime2_s2 = 4.23
J_s2 = 0.92
k_i = 0.41
k_o = 0.05
k_prime2_d2 = 0.79
k_repair = 0.08
J1 = 0.1
J2 = 0.1

# helper function
def Goldbeter_koshland(u,v,q,r):
    a = 2*u*r
    b = (v-u+v*q+u*r)**2
    c = 4*u*r*(v-u)
    d = math.sqrt(b-c)
    e = v-u+v*q+u*r+d
    return a/e

def heavisible(x):
    if x > 0:
        return 1
    else:
        return 0

# all the ODE
m = GEKKO()
m.time = np.linspace(0,40,1000)
X_var = m.Var(value=0.0)
Y_var = m.Var(value=0.19)
Z_var = m.Var(value=0.78)
P_var = m.Var(value=0.19)
k_d2 = m.Param()
k_d53 = m.Param()
m.Equation(X_var.dt()==-k_repair*heavisible(X_var))
m.Equation(k_d2==k_prime_d2*(1+X_var))
m.Equation(k_d53=k_prime_d53+k_prime2_d53*Goldbeter_koshland(Z_var,theta,J1/P_var,J2/P_var))
m.Equation(P_var.dt()==k_prime_s53+k_prime2_s53*(Y_var**4/(J_s53**4+Y_var**4))-k_prime_d53*P_var)
m.Equation(Y_var.dt()==k_prime_s2 + k_prime2_s2*(P_var**4/(J_s2**4+P_var**4))-k_i*Y_var+k_o*Z_var-k_prime_d2*Y_var)
m.Equation(Z_var.dt()==k_i*Y_var-k_o*Z_var-k_d2*Z_var)
m.options.IMODE=4
m.solve()






