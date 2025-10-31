# rational-approximation
This is a Fortran library for rational approximation of functions. The library can be built using the Fortran Package Manager (fpm). Multiple precisions are supported using the fypp preprocessing tool. Both real and complex types are supported.

## Currently supported interpolation schemes
- ```thiele_interp``` Interpolates using Thiele's interpolation formula 

## Currently supported approximation schemes

## Building the library
The library can be built with the ```fpm build``` command.

## Examples
The examples can be found in the ```example``` directory. 

To run all of them, do: ```fpm run --example```. 

To run a specific example, instead do: ```fpm run --example Name_of_example```.

The current list of examples is:
- ```thiele_sin``` Computes an interpolation of sin for real arguments using Thiele's interpolation formula (```thiele_interp```)

## Supported compilers
The library has been tested with the following compilers (on a Linux system):
- ```ifx``` version 2025.3.0

Earlier versions may work, but have not been tested at the moment. 

The library uses parameterized derived types (PDTs) with ```kind``` type parameters, which means this must be supported by the compiler you intend to build the library with.
