import numpy as np

# IC

def en(rho,P,v,gamma):
    return (P/(rho*(gamma-1))) + 0.5*v**2

def IC(x,dx,gamma,system,xd,epsilon):
    if system==0:
        rho_l, P_l, v_l = 1., 3.5, 0.
        rho_r, P_r, v_r = 0.2, 0.5, 0.
    elif system==1:
        rho_l, P_l, v_l = 1., 1e6, 0.
        rho_r, P_r, v_r = 1., 1., 0.

    rho = rho_r + 0.5*(1 - np.tanh((x-xd)/(epsilon*dx)))*(rho_l - rho_r)
    P = P_r + 0.5*(1 - np.tanh((x-xd)/(epsilon*dx)))*(P_l - P_r)
    v = v_r + 0.5*(1 - np.tanh((x-xd)/(epsilon*dx)))*(v_l - v_r)
    e = en(rho,P,v,gamma)

    cs = np.sqrt(gamma*P/rho)

    rhov = rho*v
    rhoe = rho*e

    return rho,rhov,rhoe,v,e,P,cs
