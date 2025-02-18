# Instructions to run the incremental driver
All the instructions have been tested on Linux, using gfortran as compiler.

## Compilation
In order to compile the incremental driver, you need to compile the material model and the incremental driver
as an object file. Then you need to link both object files to create the executable.


### Compile the material model
To compile the material model, e.g. the Linear Elastic material model you need to run the following command:

```bash
gfortran -c elastic.f -o elastic.o
```

### Compile the incremental driver
To compile the incremental driver, you need to run the following command:

```bash
gfortran -c incrementalDriver.f -o incrementalDriver.o
```

# Link the material model and the incremental driver
To link the material model and the incremental driver, you need to run the following command:

```bash
gfortran elastic.o incrementalDriver.o -o incrementalDriver
```

Before run you need to change the incrementalDriver permissions

```bash
chmod +x incrementalDriver
```

## Run the incremental driver with elastic umat
To run the incremental driver with the elastic material model, you need to run the following command:

```bash
./incrementalDriver test=./example/test.inp param=./example/parameters.inp ini=./example/initialconditions.inp
```

The output of the program will be stored in the file `Txt_UnComAni.out`.
In the folder [example](./example) you can find the input files for the test, the parameters and the initial conditions.

