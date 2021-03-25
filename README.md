# PHSX815 Spring 2021 Week 8

## Parameter estimation & Maximum Likelihood

This repository contains the following program(s):

- `HW8/Neyman.x` : an adaptation of the macro found at `crogan/PHSX815_Week8` to
                   run from the command line with the addition of several
                   options such as choice of distribution (between
                   pre-configured a Gaussian and a new exponential dist),
                   number of measurements per experiment, and which "slice" to
                   view. Also automatically saves the plots to picture files.
                   [C++]


### Requirements

In order to compile (by typing `make`) and run the C++ examples, you
need the ROOT package installed (for visualization):
- [ROOT](https://root.cern/) (C++)

### Usage

All of the executables can be called from the
command line with the `-h` or `--help` flag, which will print the options
