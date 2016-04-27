# Project 1

This project aims to solve Poissons equation,

\frac{du^2}{dx^2} = f(x)

for a given function f. This is achieved by employing matrix methods, and deriving an algorithm for the special case of a tridiagonal matrix.

## Usage
The executable "project1" should be run with the number of meshpoints N and a filename to save the data to along with the method chosen. The python script "plot.py" can be called with the filname used in the executable to plot the data.

Alternatively a wrapper can be used by using "launcher.py". This asks for 
the inputs and validates it before it calls the executables, and gives both the file and a plot of the solution.

## Installation
A makefile is included for easy installation. Open a terminal and type "make" to make the executables.

## Requirements
This program requires the following libraries
* armadillo
* openblas
* lapack
* superlu

If superlu is not installed the sparse methods can not be used.