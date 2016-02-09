#!/bin/bash
# for compiling and running execute the following commands

ffc -O -l dolfin Forms.ufl

cmake .

make

