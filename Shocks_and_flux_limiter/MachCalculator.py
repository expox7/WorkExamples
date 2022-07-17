import numpy as np
import matplotlib.pyplot as plt
import Limiter as lmt
import InitialConditions as IC

def M0(gamma,p0,p1):
    return np.sqrt(((gamma+1)/(2.*gamma))*((p1/p0)+(gamma-1)/(gamma+1)))

scheme = 0 # 0 -> LF, 1 -> LWR, 2 -> Variable scheme with Flux Limiter
system = 0 # 0 -> Sod, 1 -> Supernova
limiter = 'Minmod' # Minmod, Osher, VanAlbada
plot = True # True -> Plot, False -> No Plot
if system==0: # Final time
    stop = 0.1
elif system==1:
    stop = 0.001

#Parameters
alfa = 0.9 # CFL number
gamma = 5./3 # Adiabatic index
epsilon = 5 # Size of the shock
if system==0: # Position of the shock
    xd = 0.5
elif system==1:
    xd = 0.1

# Domain
xi = 0.
xf = 1.
N = 500
dx = (xf-xi)/N
x = np.arange(xi-0.5*dx, xf+1.5*dx, dx)

# IC
rho,rhov,rhoe,v,e,P,cs = IC.IC(x,dx,gamma,system,xd,epsilon)


t=0.
i=1 # Integer for percentage printing

x_shock = np.array([]) # Position and velocity of the shock
mach_an = np.array([]) # Analytical Mach numbers
time = np.array([])

while t<stop:
    
    dt = alfa*(dx/np.max(v+cs))

    rho,rhov,rhoe,v,e,P = lmt.update(rho,rhov,rhoe,v,e,P,dx,dt,gamma,scheme,limiter)

    t+=dt
    
    if t > stop/2. and limiter == 'Minmod': # At this time the shock will be fully created and propagating normally.
    
        v_dif = v[:-1] - v[1:] # Maximum value will give us the position of the shock.
        
        x_shock = np.append(x_shock,x[np.argmax(v_dif)]) # Array with the shock position at all times.
        time = np.append(time,t)
        
        mach_an = np.append(mach_an,M0(gamma,P[np.argmax(v_dif)+int(N/10)],P[np.argmax(v_dif)-int(N/10)]))
    
    if (t-dt < i*(stop/10.)) and (i*(stop/10.) < t+dt): # Prints percentage of completion. It assumes a constant dt.
        print(i*10, '%')
        i+=1

if limiter == 'Minmod':    
    print('Analytical Mach number: ', np.mean(mach_an))
    
    v_shock = np.polyfit(time,x_shock,1)[0]
    
    print('Numerical Mach number:', np.mean(v_shock/cs[-1]))

if plot:
    fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1)
    
    line1, = ax1.plot([],[])
    line2, = ax2.plot([],[])
    line3, = ax3.plot([],[])
    line4, = ax4.plot([],[])
    
    ax1.set_xlim(xi,xf)
    ax2.set_xlim(xi,xf)
    ax3.set_xlim(xi,xf)
    ax4.set_xlim(xi,xf)
    
    if system==0:
        ax1.set_ylim(0.,1.2)
        ax2.set_ylim(-3.,3.)
        ax3.set_ylim(0.,4.)
        ax4.set_ylim(2.,8.)
    elif system==1:
        ax1.set_ylim(0.,5.)
        ax2.set_ylim(0.,900.)
        ax3.set_ylim(0.,1.5e6)
        ax4.set_ylim(0.,1.6e6)
    
    ax4.set_xlabel(r'$x$')
    
    ax1.set_ylabel(r'$\rho$')
    ax2.set_ylabel(r'$v$')
    ax3.set_ylabel(r'$P$')
    ax4.set_ylabel(r'$e$')
    
    ax1.set_title(r"$t=%.6f$"%t)
    line1.set_data(x,rho)
    line2.set_data(x,v)
    line3.set_data(x,P)
    line4.set_data(x,e)

    ax1.grid()
    ax2.grid()
    ax3.grid()
    ax4.grid()
    
    plt.show()
