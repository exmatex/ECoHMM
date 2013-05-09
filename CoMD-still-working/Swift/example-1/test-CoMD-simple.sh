#!/bin/bash

stc -r ${PWD} -r ${HOME}/exm/apps/swig-data CoMD-simple.{swift,tcl}

turbine -l -n 3 CoMD-simple.tcl
