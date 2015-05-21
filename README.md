# Genetic algorithms for simulating optimal food web structures (fwga)

This README is for the R package fwga (food web genetic algorithms) (see also
Yen JDL, et al., Highly connected food webs maximize robustness but not ecosystem throughput, in review).

Copyright &copy; 2015, Jian Yen

*****

## Licence details
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*****

## Overview
fwga is a collection of R functions for simulating optimal food web (interaction network) structures according to one or more fitness functions. All functions in fwga are written and tested in R 3.1.2.

Several additional packages are required (see Installation, below) but all are available through CRAN.

Created 21 May 2015

Last updated 21 May 2015

*****

## Installation
fwga is distributed as an R package (in source form) but is not currently available through the CRAN. The current version of fwga has been tested only on OSX 10.10.3 but should be suitable for Windows and other Unix systems. A binary version of this package will be uploaded soon.

Download the file fwga_x.x.tar.gz into a local directory on your computer (replace x.x with the current version).

Because fwga is installed from source you will also need an appropriate C and C++ compiler installed. Easily installed options are [gcc](https://github.com/kennethreitz/osx-gcc-installer/) (OSX users) and [Rtools](https://github.com/stan-dev/rstan/wiki/Install-Rtools-for-Windows) (Windows users). If you're unsure of whether you need to install a C/C++ compiler, you can try installing the fwga package anyway; if you do not get any errors then no compiler is needed.

fwga imports functions from several packages (see Depends and Imports in the DESCRIPTION file) and these packages must be installed for fwga to install and load correctly. All packages are available through the CRAN and should be easy to install.

Once the appropriate packages have been installed, you simply need to place the file fwga_x.x.tar.gz in the current working directory (where x.x is replaced with the current version number) and use
```
install.packages("fwga_x.x.tar.gz", repos=NULL, type="source")
```
to install the fwga R package. This package then can be loaded in R using `library(fwga)`.

Errors during installation often are related to the installation of required packages. Restarting R and making sure all required packages can be loaded using `library(pkgName)` will identify missing or incorrectly installed packages.

## Usage
Once fwga has been installed there is one main function to use: `fwga`. This function has been set up with reliable default settings and information about its use can be found by typing `?fwga` in the R console.

The main settings to change are the number of iterations `n.iter`, the number of populations `n.pop`, the starting population `suggest`, and the desired fitness function/s `fits`. With default settings and enough iterations, the genetic algorithm should converge to optimal food-web structures. Convergence can be assessed by plotting the calculated fitness against iterations; type `?fwga` into the R console for more information. Note that the `fwga` function minimizes the chosen fitness functions, so, if maximization is required, the fitness function must be multiplied by minus one. 

Mathematical details of the genetic algorithm, along with an example of its application, are included in:
Yen JDL, et al. (in review) Highly connected food webs maximize robustness but not ecosystem throughput.


*****

## Feedback
Please send comments and bug reports to
<jdl.yen@gmail.com>
