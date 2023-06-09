# formikoj: A flexible library for data management and processing in geophysics

by
Matthias Steiner and Adrián Flores Orozco


![formikoj logo](./media/logo/formikoj_logo.svg "formikoj logo")

[![Python: 3.7](https://img.shields.io/badge/Python-3.7-blue.svg)](#)
[![License](https://img.shields.io/badge/license-MIT%20License-green)](LICENSE.md)

## Abstract

We introduce here the open-source library formikoj, which provides a flexible 
framework for managing and processing geophysical data collected in 
environmental and engineering investigations. To account for the substantial 
changes regarding the market shares of operating systems within the last two 
decades, the library is specifically implemented and tested for 
cross-plattform usage.

![Graphical abstract](./manuscript/graphical_abstract.svg)

## Structure of this repository

The source code of the formikoj library is in the `code` folder.
Data and information to reproduce the exemplary use cases presented in the 
manuscript are provided in the `examples` folder. 
The pdf and graphical abstract of the manuscript are in the `manuscript` folder.
See the `README.md` files in each directory for a full description.

## Getting the code

You can download a copy of all the files in this repository by cloning the
[git](https://git-scm.com/) repository:

    git clone https://github.com/geophilik/formikoj 

or [download a zip archive](https://github.com/geophilik/formikoj/archive/refs/heads/main.zip).

## Dependencies

You'll need a working Python environment to run the code.
The recommended way to set up your environment is through the
[Anaconda Python distribution](https://www.anaconda.com/download/) which
provides the `conda` package manager.
Anaconda can be installed in your user directory and does not interfere with
the system Python installation.
The required dependencies are specified in the file `environment.yml`.

We use `conda` virtual environments to manage the project dependencies in
isolation.
Thus, you can install our dependencies without causing conflicts with your
setup (even with different Python versions).

Open a terminal (Linux & Mac) or the Anaconda Prompt (Windows) and run the 
following command in the repository folder (where `environment.yml`
is located) to create a separate environment and install all required
dependencies in it:

    conda env create
    
If you want to use the formikoj library in an existing conda environment
you install the missing the dependencies based on `environment.yml` as
described [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html?highlight=prune#updating-an-environment)

## Installing the library

To use the formikoj library we suggest to install it in an conda 
environment. To do so you first activate the corresponding environment and then
use the Makefile provided in the `code` folder to run the setup process.

    conda activate <env_name>
    cd code
    make build

## License

All source code is made available under the MIT License. See LICENSE.md for 
the full license text.

## Credits

The structure of this reprository is based on this [template](https://github.com/pinga-lab/paper-template) 
for reproducible papers by Leonardo Uieda as well as the adaptions implemented
in [this](https://github.com/florian-wagner/four-phase-inversion) repository by 
Florian Wagner.
