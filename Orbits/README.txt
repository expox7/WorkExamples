The folder has two python files; update.py and orbits.py.

update.py is a python module with the acceleration function and two EDO solver 
function; one with the Euler method and the other with the rk4 method, prepared for
the dynamic of particles in specific. These functions take masses, positions and velocities
at certain t and output the positions and velocities at t + dt, being dt the established
timestep.

orbits.py use update.py to calculate the orbits of a certain set of particles with set
masses and initial position and velocity in an NFW potential. It can also calculate orbits
taking into consideration the mutual gravitational attraction of the particles.

The program is connected to the ICs provided in the input folder (disk10.txt, disk100.txt and disk1000.txt)
and to two other folders; Output, where it puts the snapshots if the user wants any, and AnimationOutput, 
where it puts the images and animation of the particle movement over time if the user request it.

Project_Orbits.pdf has the information about the project.
