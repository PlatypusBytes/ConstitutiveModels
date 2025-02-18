# Instructions to run the incremental driver
All the instructions have been tested on Linux, using gfortran as compiler.

## Compilation
In order to compile the incremental driver, you need to compile the material model and the incremental driver
as an object file. Then you need to link both object files to create the executable.


### Compile the material model
To compile the material model, e.g. Mohr-Coulomb material model you need to run the following command:

```bash
gfortran -c UMAT_MohrCoulomb.f -o UMAT_MohrCoulomb.o
```

### Compile the incremental driver
To compile the incremental driver, you need to run the following command:

```bash
gfortran -c incrementalDriver.f -o incrementalDriver.o
```

# Link the material model and the incremental driver
To link the material model and the incremental driver, you need to run the following command:

```bash
gfortran UMAT_MohrCoulomb.o incrementalDriver.o -o incrementalDriver
```

Before run you need to change the incrementalDriver permissions

```bash
chmod +x incrementalDriver
```

## Run the incremental driver
To run the incremental driver with the compiled material model (in this example the Mohr-Coulomb material model),
 you need to run the following command:

```bash
./incrementalDriver test=test.inp param=parameters.inp ini=initialconditions.inp
```

The output of the program will be stored in the file specified in the first line of `test.inp` file.
In the folder [example_MorhCoulomb](./example_MorhCoulomb) you can find the input files for the `test.inp`,
the `parameters.inp` and the `initialconditions.inp`.

