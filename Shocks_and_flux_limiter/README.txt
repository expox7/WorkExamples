There are four python files: InitialConditions, Limiter, ShockAnimation and MachCalculator.

InitialConditions it's a module with a function that returns the ICs when given the needed parameters.
It has the ICs for the sod shock and the supernova explosion.

Limiter it's another module with all the function needed for the limiter implementation. 
The main function is "update", which give the variables at n+1 when given the variables at n
and the necessary information (scheme, limiter...).

ShockAnimation and MachCalculator are the programs the user need to use to run the simulations.

The first one makes an animation of the evolution of the system from time zero.
The second one plots (if the user asks for it) the state of the system at a certain time and 
calculates the numerical and analytical Mach numbers of the shock.

Both of them have a list of options to manipulate at the start of the program:

scheme -> Determines if the LF, LFR or flux limiter scheme is used.
system -> Determines if the program simulates the sod shock or the supernova explosion.
limiter -> Determines which limiter is used if the flux limiter scheme is selected.

MachCalculator also has:

stop -> Determines the time where the simulation is stopped.
plot -> Determines if the result is plotted at t=stop.

And parameters:

alfa -> CFL number.
gamma -> Adiabatic index.
epsilon -> Size of the shock.
xd -> Position of the shock (by default it is separated by system, with a predefined value for each one).

Note: MachCalculator will only calculate the Mach numbers if the Minmod limiter is selected.
The other limiters have problems in the calculation due to the bad performance.