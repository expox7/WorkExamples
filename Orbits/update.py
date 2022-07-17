import numpy as np

G = 4.302e-6 # kpc (km/s)^2 Ms^-1
rho = 5932371 # Ms/kpc^3
Rs = 20 # kpc

def a_int(x,m,soft):
    x1 = x[:,0:1] #x position matrix. Axis 0: Particle. Axis 1: x position.
    x2 = x[:,1:2] #y position...
    x3 = x[:,2:3] #z position...

    dx1 = x1.T - x1 #x distances matrix. Axis 0: Particle. Axis 1: array with the distances to all the particles.
    dx2 = x2.T - x2 #y distances...
    dx3 = x3.T - x3 #z distances...

    inv_d3 = (dx1**2 + dx2**2 + dx3**2 + soft**2)**(-1.5) #1/d^3 matrix.

    #G * (dx1 * inv_d3) will give us a matrix like the ones before with
    #the accelerations/m produced to each particle from each other. If we
    #do matrix multiplication with the m array, the [G * (dx1 * inv_d3)]_ij
    #element of the matrix will multiplicate for m_j and all of the results will
    #sum up to give a_i in the correspondent coordinate. Note that the elements ii of dx1, dx2
    #and dx3 will be 0, so the acceleration produced in a particle by the gravitational
    #effect of itself is 0.
    #This calculation gives us a matrix with a form like x1, x2 and x3. We can transform it in
    #a unidimensional array with the corresponding ax, ay and az with the selection [:,None].
    ax1 = (G * (dx1 * inv_d3) @ m)[:,None] 
    ax2 = (G * (dx2 * inv_d3) @ m)[:,None]
    ax3 = (G * (dx3 * inv_d3) @ m)[:,None]

    ax = np.append(ax1,ax2,axis=1) #We append ax1, ax2 and ax3 in a way that produces a matrix like x, but with the accelerations.
    ax = np.append(ax,ax3,axis=1)
    
    return ax

def a(x,m,soft=0):
    r = np.sqrt(x[:,0]**2 + x[:,1]**2 + x[:,2]**2)
    M = 4*np.pi*rho*Rs**3*(np.log((Rs + r)/Rs) - r/(Rs + r))
    C = -G*M/r**3
    if soft==0:
        return C[:,None]*x
    else:
        return C[:,None]*x + a_int(x,m,soft)

def euler(x,v,dt,m,soft=0):

    vf = v + a(x,m,soft)*dt
    xf = x + v*dt

    return (xf,vf)

def rk4(x,v,dt,m,soft=0):

    kv_1 = a(x,m,soft)
    kr_1 = v

    kv_2 = a(x + kr_1*dt/2,m,soft)
    kr_2 = v + kv_1*dt/2

    kv_3 = a(x + kr_2*dt/2,m,soft)
    kr_3 = v + kv_2*dt/2

    kv_4 = a(x + kr_3*dt,m,soft)
    kr_4 = v + kv_3*dt
    
    vf = v + dt/6 * (kv_1+2*kv_2+2*kv_3+kv_4)
    xf = x + dt/6 * (kr_1+2*kr_2+2*kr_3+kr_4)
    
    return (xf,vf)
