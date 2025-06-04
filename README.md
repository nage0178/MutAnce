# MutAnce 
MutAnce is Bayesian inference program which estimated the time and population of origin of point mutations for structured populations using a multispecies coalescent model.
It uses the output of the program bpp an input to the program.

## Download and Compile the Program
To clone the program from GitHub

```
git clone https://github.com/nage0178/MutAnce.git
```

To compile the program, navigate to the directory and use the make command.

```
cd MutAnce/
make
```
You must have a C compiler, such as gcc or clang, installed to compile the program.
You also must have the GNU Scientific Library installed.
To use the make command, you must have the make system installed.

## Download and Compile bpp
To use the program, you must first run bpp. The example uses a development branch of bpp, so download bpp and switch branches.

```
cd ../
git clone https://github.com/bpp/bpp.git
cd bpp
git checkout mutanceDev
```

Then compile the program.

```
cd src
make
```
If you have installed the GNU scientific Library with Homebrew, use the following make command instead. 

```
make Makefile.Mac_homebrew
```

## Running the example
To run the example of MutAnce, navigate back to the MutAnce examples directory.

```
cd ../../MutAnce/examples
```

First run bpp.

```
../../bpp/src/bpp --cfile inference.ctl
```

Then run MutAnce.

```
../MutAnce example.ctl
```
## License

Copyright (C) 2025 Anna Nagel.
Parts of the code are modified from [bpp](https://github.com/bpp/bpp).

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
