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
To use the make command, you must have the make system installed.

## Download and Compile bpp
To use the program, you must first run bpp. The example uses an older version of bpp, so download bpp and switch branches.

```
cd ../
git clone https://github.com/bpp/bpp.git
cd bpp
git checkout devAnna
```

Then compile the program.

```
cd src
make
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

