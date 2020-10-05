#!/bin/bash

Qpep=
Rpep=
Qcoords=
Rcoords=
output=

/data/mg1/caix/src/biosoft/SynOrth/SynOrths -a ${Qpep} -b ${Rpep} -p ${Qcoords} -q ${Rcoords} -o ${output}
