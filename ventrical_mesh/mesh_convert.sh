#!/bin/bash

mesh_name_full=ventrical.geo
if [ -n "${1}" ]; then
    mesh_name=${1}
fi
#remove extension from file name
mesh_name_main=`echo ${mesh_name_full} | rev | cut -f 2- -d '.' | rev`

gmsh -3 -format msh2 -o ${mesh_name_main}.msh ${mesh_name}

