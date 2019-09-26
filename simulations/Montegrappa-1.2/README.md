

                     ,/k.
                    /  ih,              *MonteGrappa*
               ,-' ,  `:7b                      v1.0
             _.-/   '  /b.`.4p,
          --   ,    ,-' ^6x, `.'^=._


     WELCOME TO THE FIRST VERSION OF MONTEGRAPPA !

===========================================================================


This file will help you to compile and install Montegrappa v1.0.
The package contains three programs:

	- "montegrappa", the main Monte Carlo engine 
	- "grappino", a tool to create the input files for montegrappa.
	- "mhistogram", a data analysis tool

The code needs GSL libraries and the MPI environment to be correctly 
installed on your machine. This is a requirement to enable STEMPERING
and PARALLEL TEMPERING functions. 

Please check these dependencies before installing the software, or skip 
to the next section if you want a custom installation ( a basic version 
of Montegrappa can be compiled without GSL and MPI).

By default, GSL libraries should be installed in /opt/local/lib and their
headers installed in /opt/local/include.
If they are in different locations, you have to modify the compilation 
flags in the Makefile accordingly.

The whole distribution can be quickly compiled just typing in your terminal

	$ make all

The binary files "montegrappa", "montegrappa_mpi", "grappino" 
and "mhistogram" will be created in the ./bin directory. 




CUSTOM INSTALLATION ------------------------------------------------------


(1) MonteGrappa Single Core, no STEMPERING (no PARALLEL TEMPERING)

	$ make cleanobj
	$ make

An executable called "montegrappa" will be created in the ./bin directory.

(2) MonteGrappa Single Core, with STEPERING (no PARALLEL TEMPERING)

    - SYSTEM REQUIREMENTS: GSL libraries

	$ make cleanobj	
	$ make version=STEMPERING

As well as the (1) version, an executable called "montegrappa" will be created in the ./bin directory.


(3) MonteGrappa Multi Core, with PARALLEL TEMPERING (no STEMPERING)

    - SYSTEM REQUIREMENTS: MPI environment


	$ make cleanobj
	$ make version=MPI

An executable called "montegrappa_mpi" will be created in the ./bin directory.


(3) grappino and mhistogram --------------------------------------------------
 
The tools can be compiled alone, with the commands

      $ make grappino
      $ make mhistogram

The executables "grappino" and "mhistogram" will be created in the ./bin directory.


