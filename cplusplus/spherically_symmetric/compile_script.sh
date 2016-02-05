#!/bin/bash
# for compiling and running execute the following commands

ffc -l dolfin SphericallySymmetric.ufl

cmake .

make

