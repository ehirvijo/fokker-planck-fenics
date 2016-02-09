#!/bin/bash
# for compiling and running execute the following commands

ffc -l dolfin Forms.ufl

cmake .

make

