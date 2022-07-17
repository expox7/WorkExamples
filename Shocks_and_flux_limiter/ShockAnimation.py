import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import Limiter as lmt
import InitialConditions as IC

scheme = 0 # 0 -> LF, 1 -> LWR, 2 -> Variable scheme with Flux Limiter
system = 0 # 0 -> Sod, 1 -> Supernova
limiter = 'Minmod' # Minmod, Osher, VanAlbada

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

############## ANIMATION

if system==0:
    dt_show = 1e-3
elif system==1:
    dt_show = 1e-6

t=0.
next_show = dt_show + t

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

ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()

def animate(i):
    global t, next_show, dt, rho, rhov, rhoe, v, e, P
    
    while t<next_show:
        dt = alfa*(dx/np.max(v+cs))

        rho,rhov,rhoe,v,e,P = lmt.update(rho,rhov,rhoe,v,e,P,dx,dt,gamma,scheme,limiter)
        
        t+=dt

    next_show = t+dt_show

    ax1.set_title(r"$t=%.5f$"%t)
    line1.set_data(x,rho)
    line2.set_data(x,v)
    line3.set_data(x,P)
    line4.set_data(x,e)

# This function takes as arguments:
#   * fig : the figure we will use for plotting purposes
#   * animate : the function that will be called each frame
#   * interval : time in ms between frames
#   * blit : clear the figure each frame
#   * fargs : extra arguments to 'animate'
#
# Please read carefully the documentation
ani = animation.FuncAnimation(fig, animate, interval=100, blit=False)
plt.show()
