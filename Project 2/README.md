#Project 2

This project solves the Schrödinger equation by finding the energy levels.

The problem is reduced to the matrix equation
Ax = \lambda x
which is solved as an eigenvalue-problem with Jacobi rotations.

##Installation
Open the terminal in this folder and type

```
make
```

##Usage
In the terminal type

```
./project2
```

This gives the input parameters and executes the unit tests. These parameters are


N, the size of the system. This should be a number less than 300 for a reasonable computer.

rho_inf, the maximum value for the axis. If the resulting eigenvector "collides" with this axis, or the eigenvalues are incorrect, this must be chosen larger.

eig_vec, a flag which determines if the eigenvalues should be calculated or not,
0 for false.

omega_r, the parameter for the coupling between two electrons in the harmonic oscillator. If omitted, the single electron in the harmonic oscillator is calculated.

The program returns a "test.txt" containing the three lowest energies and the corresponding eigenvectors.


To plot the eigenvectors a python script is included. Usage

```
python3 plot.py
```

when a textfile has been generated with the above program. These can be chained to give the plot straight away, for example;

```
./project2 100 6 1 && python3 plot.py
```
Which will plot the three lowest energies and the associated eigenvectors (squared).