#!/bin/bash
# for compiling and running execute the following commands

ffc -O -l dolfin Forms.ufl

cmake .

make

if [ $? -eq 0 ]; then
    echo "Generating mesh"
    ./generate_fenics_mesh.sh
else
    echo "C++ compiling errors. Verify your code."
fi


