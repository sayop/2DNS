Code development
================

The present project is aimed to develop a computer program for solving 2-D unsteady Navier-Stokes equations for a supersonic problem. Hereafter, the program developed here in this project is called '2DNS'.

2DNS Code summary
-----------------
The source code contains two directories, 'io', and 'main', for input/output related sources and main solver routines, respectively. 'CMakeLists.txt' file is also included for cmake compiling.

::

   $ cd 2DNS/CODEdev/src/
   $ ls
   $ CMakeLists.txt  io  main

The **io** folder has **io.F90** and **ReadGrid.F90** files which contains subroutines for reading input/output data and grid info. It also includes **input** directory which contains a default **input.dat** file.

The **main** folder is only used for calculating essential subroutines required to solve the '2DNS' equation by using 'AUSMPW+' scheme for solving inviscid flux reconstruction and independent viscous flux calculator. The main routine is run by **main.F90** which calls important subroutines from **main** folder itself and **io** folder when needed.

