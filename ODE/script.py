import numpy as np
from scipy.integrate import odeint
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from gekko import GEKKO

'''
Figure1A linear signal response
dR/dt = k0 + k1*S - k2*R
'''

def model(R,t,S):
    k0 = 0.01
    k1 = 1
    k2 = 5
    dRdt = k0 + k1*S - k2*R
    return dRdt

R0 = [0,0.3,0.5]
t = np.linspace(0,1,10)
S = 1
result = odeint(model,R0,t,args=(S,))

# solution curve
fig,ax = plt.subplots()
ax.plot(t,result[:,0],label='R0=0')
ax.plot(t,result[:,1],label='R0=0.3')
ax.plot(t,result[:,2],label='R0=0.5')
ax.legend()
ax.set_xlabel('t')
ax.set_ylabel('R')
ax.axhline(y=0.202,xmin=0,xmax=1,linestyle='--',c='k')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1A_solution_curve.pdf',bbox_inches='tight')
plt.close()

# rate curve
k0 = 0.01
k1 = 1
k2 = 5
fig,ax = plt.subplots()
S_options = [1,2,3]
for S in S_options:
    R = np.linspace(0,1,10)
    removal_rate = k2 * R
    production_rate = [k1 * S] * len(R)
    ax.plot(R,removal_rate,linestyle='-',c='k')
    ax.plot(R,production_rate,linestyle='--',c='k')
ax.set_xlim(0,1)
ax.set_ylim(0,6)
ax.set_xlabel('R')
ax.set_ylabel('Rate')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1A_rate_curve.pdf',bbox_inches='tight')
plt.close()

# signal-response curve
k0 = 0.01
k1 = 1
k2 = 5
S = np.linspace(0,3,7)
R_ss = (k0 + k1*S) / k2
fig,ax = plt.subplots()
ax.plot(S,R_ss,linestyle='-',c='k')
ax.set_xlim(0,3)
ax.set_ylim(0,0.7)
ax.set_xlabel('Signal(S)')
ax.set_ylabel('Response(R)')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1A_SR_curve.pdf',bbox_inches='tight')
plt.close()

'''
Figure1B, hyperbolic response
dRp/dt = k1*S(Rt-Rp) - k2*Rp
'''

def model(y,t,S):
    k1 = 1
    k2 = 1
    Rt = 1
    dydt = k1*S*(Rt-y) - k2*y
    return dydt

S = 1
Rp0 = [0,0.5,1]
t = np.linspace(0,1,10)
result = odeint(model,Rp0,t,args=(S,))


# solution curve
fig,ax = plt.subplots()
ax.plot(t,result[:,0],label='R0=0')
ax.plot(t,result[:,1],label='R0=0.5')
ax.plot(t,result[:,2],label='R0=1')
ax.legend()
ax.set_xlabel('t')
ax.set_ylabel('Rp')
ax.axhline(y=0.5,xmin=0,xmax=1,linestyle='--',c='k')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1B_solution_curve.pdf',bbox_inches='tight')
plt.close()




# signal-response curve
k1 = 1
k2 = 1
Rt = 1
S = np.linspace(0,10,100)
Rp_ss = (S * Rt) / (k2/k1 + S)
fig,ax = plt.subplots()
ax.plot(S,Rp_ss,linestyle='-',c='k')
ax.set_xlim(0,10)
ax.set_ylim(0,1.1)
ax.set_xlabel('Signal(S)')
ax.set_ylabel('Response(R)')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1B_SR_curve.pdf',bbox_inches='tight')
plt.close()

'''
sigmoidal curve
dRp/dt = (k1*S*(Rt-Rp)/(km1+Rt-Rp)) - k2*Rp/(km2+Rp)
'''

def model(Rp,t,S):
    k1 = 1
    k2 = 1
    Rt = 1
    km1 = 0.05
    km2 = 0.05
    dRpdt = (k1*S*(Rt-Rp)/(km1+Rt-Rp)) - k2*Rp/(km2+Rp)
    return dRpdt

S = 1
Rp0 = [0,0.3,1]
t = np.linspace(0,20,200)
result = odeint(model,Rp0,t,args=(S,))

# solution curve
fig,ax = plt.subplots()
ax.plot(t,result[:,0],label='R0=0')
ax.plot(t,result[:,1],label='R0=0.3')
ax.plot(t,result[:,2],label='R0=1')
ax.legend()
ax.set_xlabel('t')
ax.set_ylabel('Rp')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1C_solution_curve.pdf',bbox_inches='tight')
plt.close()


# rate curve
k1 = 1
k2 = 1
Rt = 1
km1 = 0.05
km2 = 0.05
Rp = np.linspace(0,1,100)
fig,ax = plt.subplots()
for S in [0.25,0.5,1,1.5,2]:
    removal_rate = k2*Rp/(km2+Rp)
    production_rate = k1*S*(Rt-Rp)/(km1+Rt-Rp)
    ax.plot(Rp,removal_rate,linestyle='-',c='k')
    ax.plot(Rp,production_rate,linestyle='--',c='k')
ax.set_xlim(0,1)
ax.set_xlabel('Rp')
ax.set_ylim(0,2)
ax.set_ylabel('Rate')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1C_rate_curve.pdf',bbox_inches='tight')
plt.close()

# signal-response curve
S_all = np.linspace(0,3,100)
def equation(Rp,S):
    k1 = 1
    k2 = 1
    Rt = 1
    km1 = 0.05
    km2 = 0.05
    return k1*S*(Rt-Rp)/(km1+Rt-Rp) - k2*Rp/(km2+Rp)

from scipy.optimize import fsolve
store = []
for S in S_all:
    Rp_ss = fsolve(equation,[1],args=(S,))[0]
    store.append(Rp_ss)

fig,ax = plt.subplots()
ax.plot(S_all,store,c='k')
ax.set_xlim(0,3)
ax.set_xlabel('Signal(S)')
ax.set_ylim(0,1.1)
ax.set_ylabel('Response(R_ss)')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1C_SR_curve.pdf',bbox_inches='tight')


'''
Figure1d perfectly adaptation to signal
dR/dt = k1*S - k2*X*R
dX/dt = k3*S - k4*X
'''

# solution curve
from gekko import GEKKO
k1,k2,k3,k4 = 2,2,1,1
m = GEKKO()
m.time = np.linspace(0,20,201)
S = np.array([0] * 40 + [1] * 40 + [2] * 40 + [3] * 40 + [4] * 41)
S_param = m.Param(value=S)
X_var = m.Var(value=0)
R_var = m.Var(value=0)
m.Equation(X_var.dt()==k3*S_param-k4*X_var)
m.Equation(R_var.dt()==k1*S_param-k2*X_var*R_var)
m.options.IMODE = 4
m.solve()

fig,ax1 = plt.subplots()
ax1.plot(m.time,R_var,c='k')
ax1.plot(m.time,X_var,c='g')
ax1.set_xlim(0,20)
ax1.set_xlabel('Time')
ax2 = ax1.twinx()
ax2.plot(m.time,S_param,c='r')
import matplotlib.lines as mlines
ax1.legend(handles=[mlines.Line2D([],[],c=i,linestyle='-') for i in 'kgr'],labels=['R','X','S'])
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1d_solution_curve.pdf',bbox_inches='tight')
plt.close()

# rate curve
S_options = [1,2,3]
fig,ax = plt.subplots()
for S in S_options:
    R = np.linspace(0,2,100)
    X_ss = k3*S/k4
    removal_rate = k2 * X_ss * R
    production_rate = [k1*S] * len(R)
    tmp_dict = {1:'r',2:'purple',3:'g'}
    ax.plot(R,removal_rate,c=tmp_dict[S],linestyle='-')
    ax.plot(R,production_rate,c=tmp_dict[S],linestyle='--')
ax.set_xlim(0,2)
ax.set_ylim(0,8)
ax.set_xlabel('R')
ax.set_ylabel('Rate')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1d_rate_curve.pdf',bbox_inches='tight')
plt.close()

'''
Figure1e mutual activation, one-way switch
'''
k0,k1,k2,k3,k4,J3,J4 = 0.4,0.01,1,1,0.2,0.05,0.05
def Goldbeter_Koshland(v1,v2,J1,J2):
    '''
    will return the equalibrium concertration of either phosphorated or unphosphorated form of E,
    v1,J1: incoming
    v2,J2: outcoming
    '''
    B = v2-v1+J1*v2+J2*v1
    equilibrium = 2*v1*J2/(B+np.sqrt(B**2-4*(v2-v1)*v1*J2))
    return equilibrium

def model(R,t,S):
    EP = Goldbeter_Koshland(k3*R,k4,J3,J4)
    dRdt = k0*EP + k1*S - k2*R
    return dRdt

t = np.linspace(0,10,100)
result = odeint(model,y0=1,t=t,args=(0,))

# solution curve, they don't ask me to do that, skip for now

# rate curve
R = np.linspace(0,0.7,100)
S_options = [0,8,16]
fig,ax = plt.subplots()
for S in S_options:
    removal_rate = k2*R
    EP = Goldbeter_Koshland(k3*R,k4,J3,J4)
    production_rate = k0*EP+k1*S
    ax.plot(R,removal_rate,c='k',linestyle='-')
    ax.plot(R,production_rate,c='k',linestyle='--')
ax.set_xlim(0,0.7)
ax.set_xlabel('R')
ax.set_ylim(0,0.6)
ax.set_ylabel('Rate')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1e_rate_curve.pdf',bbox_inches='tight')
plt.close()

# signal-response curve
S_options = np.linspace(0,15,100)
def equation(R,S):
    EP = Goldbeter_Koshland(k3*R,k4,J3,J4)
    return k0*EP + k1*S - k2*R

# it's necessary to first check the vector field
fig,ax = plt.subplots()
S = 10
test = [equation(R,S) for R in np.linspace(0,1,100)]
ax.plot(np.linspace(0,1,100),test)
# you can find two stable steady state and an unstable steady state

fig,ax = plt.subplots()
store_uplimb,store_downlimb = [],[]
for S in S_options:
    R_ss_uplimb = fsolve(func=equation,x0=[1],args=(S,))[0]
    R_ss_downlimb = fsolve(func=equation,x0=[0],args=(S,))[0]
    store_uplimb.append(R_ss_uplimb)
    store_downlimb.append(R_ss_downlimb)
ax.plot(S_options,store_uplimb,c='k')
ax.plot(S_options[0:67],store_downlimb[0:67],c='k')
ax.plot(np.linspace(0,10,5),[0.2,0.18,0.17,0.16,0.14],c='k',linestyle='--')
ax.set_xlim(0,15)
ax.set_ylim(0,0.7)
ax.set_xlabel('Signal(S)')
ax.set_ylabel('Response(R)')
ax.text(x=7.5,y=0.65,s='Mutual activation',horizontalalignment='center',verticalalignment='center')
plt.savefig('/Users/ligk2e/Desktop/ODE/Figure1e_SR_curve.pdf',bbox_inches='tight')


'''
Figure1f mutual inhibition
'''
k0=0
k1=0.05
k2=0.1
k2_prime=0.5
k3=1
k4=0.2
J3=0.05
J4=0.05

# rate curve
def Goldbeter_Koshland(v1,v2,J1,J2):
    '''
    will return the equalibrium concertration of either phosphorated or unphosphorated form of E,
    v1,J1: incoming
    v2,J2: outcoming
    '''
    B = v2-v1+J1*v2+J2*v1
    equilibrium = 2*v1*J2/(B+np.sqrt(B**2-4*(v2-v1)*v1*J2))
    return equilibrium
R = np.linspace(0,1.5,100)
S_options = [0.6,1.2,1.8]
ER = Goldbeter_Koshland(k4, k3 * R, J3, J4)
removal_rate = k2 * R + k2_prime * ER * R
fig,ax = plt.subplots()
ax.plot(R,removal_rate,color='k',linestyle='-')
for S in S_options:
    production_rate = [k0+k1*S] * 100
    ax.plot(R,production_rate,color='k',linestyle='--')
ax.set_xlim(0,1.5)
ax.set_ylim(0,0.1)
ax.set_xlabel('R')
ax.set_ylabel('Rate(dR/dt)')
ax.text(x=1.3,y=0.035,s='0.6',horizontalalignment='center',verticalalignment='center')
ax.text(x=1.3,y=0.065,s='1.2',horizontalalignment='center',verticalalignment='center')
ax.text(x=1.3,y=0.095,s='1.8',horizontalalignment='center',verticalalignment='center')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1f_rate_curve.pdf',bbox_inches='tight')
plt.close()

# SR curve
S_options = np.linspace(0,2,100)
def equation(R,S):
    ER = Goldbeter_Koshland(k4, k3 * R, J3, J4)
    return k0+k1*S-k2*R-k2_prime*ER*R
# peek this equation, traverse starts from 0 or from 1 will end up with two stable solution, one unstable in between
S = 1
R_options = np.linspace(0,1,100)
result = [equation(R,S) for R in R_options]
fig,ax=plt.subplots()
ax.plot(R,result)
# start to plot
fig,ax=plt.subplots()
uplimb,downlimb = [],[]
for S in S_options:
    uplimb.append(fsolve(func=equation,x0=[1],args=(S,))[0])
    downlimb.append(fsolve(func=equation,x0=[0],args=(S,))[0])
ax.plot(S_options[round(100*0.84/2):],uplimb[round(100*0.84/2):],c='k')
ax.plot(S_options[:round(100*1.7/2)],downlimb[:round(100*1.7/2)],c='k')
ax.plot(np.linspace(0.842,1.7,5),[0.285,0.25,0.23,0.20,0.17],c='k',linestyle='--')
ax.set_xlim(0,2)
ax.set_ylim(0,1)
ax.set_xlabel('Signal(S)')
ax.set_ylabel('Response(R)')
ax.set_title('Mutual inhibition')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1f_SR_curve.pdf',bbox_inches='tight')
plt.close()


'''
figure1g homestatis, negative feedback
'''
k0=1
k2=1
k3=0.5
k4=1
J3=0.01
J4=0.01

# rate curve
R_options = np.linspace(0,1,100)
S_options = [0.5,1,1.5]
ER = Goldbeter_Koshland(k3,k4*R,J3,J4)  # now k3 is the rate (E -> ER), damn
production_rate = k0*ER
fig,ax=plt.subplots()
ax.plot(R_options,production_rate,c='k',linestyle='--')
for S in S_options:
    removal_rate = k2*S*R
    ax.plot(R_options,removal_rate,c='k',linestyle='-')
ax.set_xlim(0,1)
ax.set_ylim(0,1.2)
ax.set_xlabel('R')
ax.set_ylabel('Rate(dR/dt)')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1g_rate_curve.pdf',bbox_inches='tight')
plt.close()

# SR curve
S_options=np.linspace(0.01,2,100)
def equation(R,S):
    ER=Goldbeter_Koshland(k3,k4*R,J3,J4)
    return k0*ER - k2*S*R
# for some reason, the fsolve can not approximate the root of above equation if you set x=1, (asymptotic issue)
points = []
for S in S_options:
    points.append(fsolve(func=equation,x0=[1],args=(S,))[0])
fig,ax = plt.subplots()
ax.plot(S_options,points,c='k')
ax.set_xlim(0,2)
ax.set_ylim(0,1)
ax.set_xlabel('Signal(S)')
ax.set_ylabel('Response(R)')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1g_SR_curve.pdf',bbox_inches='tight')
plt.close()

'''
figure2a
'''
m = GEKKO()
m.time = np.linspace(0,50,100)
#S = m.Param(np.linspace(2,5,100))
S = 4

k0 = 0
k1 = 1
k2 = 0.01
k2_prime = 10
k3 = 0.1
k4 = 0.2
k5 = 0.1
k6 = 0.05
Yt=1
Rt=1
km3=0.01
km4=0.01
km5=0.01
km6=0.01

X = m.Var(value=0)
Yp = m.Var(value=0)
Rp = m.Var(value=0)

m.Equation(X.dt() == k0+k1*S-k2*X+k2_prime*Rp*X)
m.Equation(Yp.dt() == (k3*X*(Yt-Yp))/(km3+Yt-Yp) - (k4*Yp)/(km4+Yp))
m.Equation(Rp.dt() == (k5*Yp*(Rt-Rp))/(km5+Rt-Rp) - (k6*Rp)/(km6+Rp))

m.options.IMODE=4
m.solve(disp=True)
















