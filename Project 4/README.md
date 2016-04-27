#Project 4
This project aims to find the critical temperature for the two-dimensional Ising model. The energy capacity and magnetic susceptibility is determined for a small range of values around this critical temperature.


##Installation
Open the terminal in this folder and type

```
make
```

This makes two executables, 'project3' and 'testit'. If documents and a program to animate the system is wanted run

```
make all
```

which now calls Doxygen to generate the documentation, and 

##Usage
In the terminal run 'project4', this should give the input parameters required to run the program. To run the unit test the executable 'testit' must be run,

For a visual look at the evolution of the system run 'animation', which can be specified with the options beta and L by the command line.

To reproduce the data some data has been run with a constant seed. Look under the 
Reproduction folder to see how to recreate these data-series.

