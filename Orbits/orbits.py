import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import time
import update
import imageio

#INTERNAL UNITS
##V = 1e4 # cm/s -> km/s
##M = 1.989e33 # g -> Ms
##L = 3.085678e21 # cm -> kpc

#############################   CHANGE HERE WHATEVER YOU WANT 

T = 2. # Simulated time (1 u = 1.8*10^9 years).
N = 1000 # Timesteps.

soft = 1 #Softening lenght (kpc). 0 for not self-gravitation.
euler = False #Method of ODE calculations. True for euler, False for RK4.
particles = 10 #Number of particles. Only 1, 10, 100 or 1000.

n = 0 #Number of stars plotted. 0 for all of them.

plot = True #If True, it plots the final orbits.

anim = True #If True, it makes an animation. Output in AnimationOutput folder.
if anim:
    filenames = []
    images = []
p = 20 #For animations. How many points 'to the past' plotted per star since the actual t point. The first is the actual position. 0 for all of them.

snapshots = 10 #Number of equally time-spaced snapshots. Same format as ICs. Output in Output folder. It has to be a divisor of N.

#############################

dt = T/N #Timestep length.

if particles!=1:
    if particles==10:
        f = open("inputs/disk10.txt", "r")
    elif particles==100:
        f = open("inputs/disk100.txt", "r")
    elif particles==1000:
        f = open("inputs/disk1000.txt", "r")

    x = [] #Position matrix. Axis 0: Particle. Axis 1: Position vector for each t. Axis 2: Coordinates of the position vector (x,y,z).
    v = [] #Velocity matrix. Axis 0: Particle. Axis 1: Velocity vector for each t. Axis 2: Coordinates of the velocity vector (vx,vy,vz).
    m = [] #Mass array. Axis 0: Particle.
    for line in f:
        line = line.split()
        x.append([[float(line[1]),float(line[2]),float(line[3])]])
        v.append([[float(line[4]),float(line[5]),float(line[6])]])
        m.append(float(line[0]))
    x = np.array(x)
    v = np.array(v)
    m = np.array(m)

    f.close()
else:
    x = np.array([[[8,0,0]]])
    v = np.array([[[0,127,0]]])
    m = np.array([1]) #For 1 particle, we only want to simulate the sun, not the galaxy.

print('Number of particles: ' + str(particles))
if soft==0:
    print('Without self-gravitation.')
else:
    print('Softening lenght: ' + str(soft) + ' kpc')
print('Total time: ' + str(T*1.8) + ' Gyr')
print()

start = time.time()

frame = 0
for i in range(N):
    if euler: #Here we update x and v matrix including in the Axis 1 of the position and velocity matrices the vectors corresponding to the next timestep.
        pos, vel = update.euler(x[:,-1],v[:,-1],dt,m,soft)
    else:
        pos, vel = update.rk4(x[:,-1],v[:,-1],dt,m,soft)
    x = np.append(x,pos[:,None],axis=1)
    v = np.append(v,vel[:,None],axis=1)

    if (anim) and (i>p) and (i*100./N == int(i*100./N)): #We use around 100 images to animate the evolution.

        frame += 1

        pb = int(p/2)

        #We divide the three-axis position matrix in three dual-axis matrices for each spatial coordinate and we take only the wanted points.
        #The "actual" point will be black for 1 particle simulations and density-dependent for >1 particles. The rest will be grey.
        #We flatten the resulting dual-axis matrix because we don't need to differenciate particles. We will plot all the points.
        posx_b, posy_b, posz_b = x[:,:,0][-n:,-1:].flatten(), x[:,:,1][-n:,-1:].flatten(), x[:,:,2][-n:,-1:].flatten()
        posx_g, posy_g, posz_g = x[:,:,0][-n:,-p:-1].flatten(), x[:,:,1][-n:,-p:-1].flatten(), x[:,:,2][-n:,-p:-1].flatten()

        if particles > 1:
            xyz = np.vstack([posx_b,posy_b,posz_b])
            density = stats.gaussian_kde(xyz)(xyz)
            density = density/np.amax(density)

            density = np.log10(density)

        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        ax.scatter(posx_g,posy_g,posz_g,s = 1, c='grey',marker='.')
        if particles > 1:
            ax.scatter(posx_b,posy_b,posz_b,c=density,marker='.')
        else:
            ax.scatter(posx_b,posy_b,posz_b,s = 1, c='black',marker='o')

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        ax.set_title('{:.1f}'.format(i*dt*1.8) + ' Gyr')

        if particles==1: #These numbers adapts well to each system.
            dist = 10
        elif particles==10:
            dist = 30
        elif particles==100:
            dist = 50
        elif particles==1000:
            dist = 60

        ax.set_xlim(-dist,dist)
        ax.set_ylim(-dist,dist)
        ax.set_zlim(-dist,dist)

        plt.grid()

        plt.savefig('AnimationOutput/' + str(frame) + '.png',dpi=300)
        filenames.append('AnimationOutput/' + str(frame) + '.png')          

        plt.close()

    if (snapshots != 0) and (i*snapshots/N == int(i*snapshots/N)):
        f = open('Output/snapshot' + str(int(i*snapshots/N)+1) + '.txt', 'w')
        for part in range(particles):
            f.write(str(m[part]) + ' ' + str(x[part,-1,0]) + ' ' + str(x[part,-1,1]) + ' ' + str(x[part,-1,2]) + ' ' + str(v[part,-1,0]) + ' ' + str(v[part,-1,1]) + ' ' + str(v[part,-1,2]) + '\n')
        f.close()
    
    if i*20./N == int(i*20./N): #In 5% intervals, we plot the percentage of completeness.
        print(5*i*20./N, '%')

end = time.time()

print('Computation time:', end-start)

if anim:
    for filename in filenames:
        images.append(imageio.imread(filename))
    imageio.mimsave('AnimationOutput/evolution.gif', images)

if plot:
    posx, posy, posz = x[:,:,0][-n:,:].flatten(), x[:,:,1][-n:,:].flatten(), x[:,:,2][-n:,:].flatten()
    
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.scatter(posx,posy,posz,s = 1, c='black',marker='.') #For the simple plots I prefer not to slow down the program computing the point densities.

    ax.set_xlabel('X (kpc)')
    ax.set_ylabel('Y (kpc)')
    ax.set_zlabel('Z (kpc)')

    #A cubic bounding box to simulate equal aspect ratio.
    max_range = np.array([posx.max()-posx.min(), posy.max()-posy.min(), posz.max()-posz.min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(posx.max()+posx.min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(posy.max()+posy.min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(posz.max()+posz.min())
    for xb, yb, zb in zip(Xb, Yb, Zb):
       ax.plot([xb], [yb], [zb], 'w')
    plt.grid()

    plt.figure()
    plt.plot(posx,posy,'.',color='black',markersize=1)
    plt.xlabel('X (kpc)')
    plt.ylabel('Y (kpc)')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid()

    plt.figure()
    plt.plot(posy,posz,'.',color='black',markersize=1)
    plt.xlabel('Y (kpc)')
    plt.ylabel('Z (kpc)')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid()

    plt.figure()
    plt.plot(posx,posz,'.',color='black',markersize=1)
    plt.xlabel('X (kpc)')
    plt.ylabel('Z (kpc)')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid()
    
    plt.show()
