#!/bin/bash


function main () {

    ${PYTHON_ENV}/bin/python3.8 /web.py 

}


export -f main
export TMPDIR=/tmp

export PYTHON_ENV=/opt/conda/envs/python_env



main 






