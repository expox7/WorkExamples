import numpy as np

def r(u,scheme,maxvalue=1e6):
    if scheme==0:
        return np.full(np.size(u)-2,0)
    elif scheme==1:
        return np.full(np.size(u)-2,1)
    else:
        num = u[1:-1] - u[:-2]
        den = u[2:] - u[1:-1]
        f = np.nan_to_num(num/den, nan=1.)
        f[f>maxvalue] = maxvalue
        f[f<-maxvalue] = -maxvalue
        return f

def phi(r,limiter):
    f = np.zeros(np.size(r))
    
    if limiter=='Minmod':
        f = np.maximum(f,np.minimum(np.full(np.size(r),1),r))
        
    elif limiter=='Osher':
        f = np.maximum(f,np.minimum(np.full(np.size(r),2),r))
        
    elif limiter=='VanAlbada':
        f = (r**2 + r)/(r**2 + 1)
        
    return f

def F(r,fpl,fnl,fph,fnh,limiter):
    Fp = (1-phi(r,limiter))*fpl + phi(r,limiter)*fph
    Fn = (1-phi(r,limiter))*fnl + phi(r,limiter)*fnh
    return Fp, Fn

def update(rho,rhov,rhoe,v,e,P,dx,dt,gamma,scheme,limiter):
    
    # Relative gradients
    r_rho = r(rho,scheme)
    r_rhov = r(rhov,scheme)
    r_rhoe = r(rhoe,scheme)

    ####
    
    f_rho_pl = 0.5*rhov[2:] - 0.5*dx/dt*(rho[2:] - rho[1:-1])
    f_rho_nl = 0.5*rhov[:-2] + 0.5*dx/dt*(rho[:-2] - rho[1:-1])

    f_rhov_pl = 0.5*(rhov[2:]*v[2:] + P[2:]) - 0.5*dx/dt*(rhov[2:] - rhov[1:-1])
    f_rhov_nl = 0.5*(rhov[:-2]*v[:-2] + P[:-2]) + 0.5*dx/dt*(rhov[:-2] - rhov[1:-1])

    f_rhoe_pl = 0.5*(rhoe[2:] + P[2:])*v[2:] - 0.5*dx/dt*(rhoe[2:] - rhoe[1:-1])
    f_rhoe_nl = 0.5*(rhoe[:-2] + P[:-2])*v[:-2] + 0.5*dx/dt*(rhoe[:-2] - rhoe[1:-1])

    ####
            
    rho_m = 0.5*(rho[:-1] + rho[1:]) - 0.5*dt/dx*(rhov[1:] - rhov[:-1])
    rhov_m = 0.5*(rhov[:-1] + rhov[1:]) - 0.5*dt/dx*((rhov[1:]*v[1:] + P[1:]) - (rhov[:-1]*v[:-1] + P[:-1]))
    rhoe_m = 0.5*(rhoe[:-1] + rhoe[1:]) - 0.5*dt/dx*((rhoe[1:] + P[1:])*v[1:] - (rhoe[:-1] + P[:-1])*v[:-1])

    v_m = rhov_m/rho_m
    e_m = rhoe_m/rho_m
    P_m = (e_m - 0.5*v_m**2)*rho_m*(gamma-1.)

    f_rho_ph = rhov_m[1:]
    f_rho_nh = rhov_m[:-1]

    f_rhov_ph = rhov_m[1:]*v_m[1:] + P_m[1:]
    f_rhov_nh = rhov_m[:-1]*v_m[:-1] + P_m[:-1]

    f_rhoe_ph = (rhoe_m[1:] + P_m[1:])*v_m[1:]
    f_rhoe_nh = (rhoe_m[:-1] + P_m[:-1])*v_m[:-1]

    ####

    F_rho_p, F_rho_n = F(r_rho,f_rho_pl,f_rho_nl,f_rho_ph,f_rho_nh,limiter)
    F_rhov_p, F_rhov_n = F(r_rhov,f_rhov_pl,f_rhov_nl,f_rhov_ph,f_rhov_nh,limiter)
    F_rhoe_p, F_rhoe_n = F(r_rhoe,f_rhoe_pl,f_rhoe_nl,f_rhoe_ph,f_rhoe_nh,limiter)

    ####

    rho[1:-1] = rho[1:-1] - dt/dx*(F_rho_p - F_rho_n)
    rhov[1:-1] = rhov[1:-1] - dt/dx*(F_rhov_p - F_rhov_n)
    rhoe[1:-1] = rhoe[1:-1] - dt/dx*(F_rhoe_p - F_rhoe_n)

            
    # Boundary conditions 
    rho[0] = rho[1]
    rho[-1] = rho[-2]
    rhov[0] = rhov[1]
    rhov[-1] = rhov[-2]
    rhoe[0] = rhoe[1]
    rhoe[-1] = rhoe[-2]


    # Primitive variables
    v = rhov/rho
    e = rhoe/rho
    P = (e - 0.5*v**2)*rho*(gamma-1.)

    return rho,rhov,rhoe,v,e,P
    
