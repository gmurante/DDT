# DDT
Domain Decomposition and Treebuild/walk from Gadget3

GM 12/5/2017

In DomainDecomp-and-Tree you find all the needed C and header
file, plus a Config.sh and a Makefile with pico and galileo
configured (with my account - change what you need)

In Data you find the parameter file and two GA file (initial conditions
and final snapshot).

DDT makes and initial domain decomposition and tree build.
It writes the files domaindecomp%03d.dat for each MPI task used,
that contain the domain decomposition.

Then, it looks for the neighbours of a given number n of random,
type 1 or 2 particles. Inside run.c you can change the value of n

The file look_around.c contains the driver and the function do
to a neighbours research for ONE particle in a cube of a given size.
There is a flag, NO_OUTPUT, that controls the ASCII output of
the neighbours of each particle. It is currently defined - so no 
output is done, otherwise you will have n*nTask files named 
part#######_task##.dat with the ID and task owner.

Finally, it outputs the time measured inside the code.

I took a look at the domain decomposition obtained and neightbours
found, and they seem to be completely reasonable.


<a href="https://zenodo.org/badge/latestdoi/91093112"><img src="https://zenodo.org/badge/91093112.svg" alt="DOI"></a>
