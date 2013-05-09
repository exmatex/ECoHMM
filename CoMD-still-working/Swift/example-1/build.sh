#!/bin/bash

# Build the CoMD leaf package

# This could be a Makefile but I think it is better
# to use bash as a reference example. -Justin


CoMD_HOME=$( cd $( dirname $0 ) ; cd ../.. ; /bin/pwd )
# CoMD_HOME=${HOME}/collab/CoMD

LEAF_PKG="comdswift"
LEAF_I=${LEAF_PKG}.i
LEAF_C=${LEAF_PKG}_wrap.c
LEAF_O=${LEAF_C%.c}.o
LEAF_SO="lib${LEAF_PKG}.so"
LEAF_TCL="${LEAF_PKG}.tcl"
LEAF_VERSION="0.0"

USER_C=${LEAF_PKG}.c
USER_O=${LEAF_PKG}.o

check()
{
  CODE=${?}
  if [[ ${CODE} != 0 ]]
  then
    MSG=$1
    echo ${MSG}
    exit ${CODE}
  fi
  return 0
}

# set -x

[[ -d ${CoMD_HOME} ]]
check "CoMD_HOME=${CoMD_HOME} ... not found!"

TCLSH=$( which tclsh )
check "Could not find tclsh in PATH!"

TCL_HOME=$( cd $( dirname ${TCLSH} )/.. ; /bin/pwd )
check "Could not find Tcl installation!"

echo "using Tcl in ${TCL_HOME}"

TCL_CONFIG=${TCL_HOME}/lib/tclConfig.sh

[[ -f ${TCL_CONFIG} ]]
check "Could not read tclConfig.sh!"

# This loads many Tcl configuration variables
source ${TCL_CONFIG}
check "tclConfig.sh failed!"

# Create the Tcl extension
echo "running SWIG..."
swig -tcl -module ${LEAF_PKG} ${LEAF_I}
check

CFLAGS="-Wall -fPIC -g -std=gnu99"
CFLAGS+=" -I ${CoMD_HOME}/src-flat/include"
CFLAGS+=" -I ${HOME}/exm/apps/swig-data"

# Run STC
stc -u -r ${PWD} CoMD-simple.swift CoMD-simple.tcl

# Compile the user code
echo "compile user code..."
gcc -c ${CFLAGS} ${USER_C}
check

# Compile the Tcl extension
echo "compile extension..."
gcc -c ${CFLAGS} ${TCL_INCLUDE_SPEC} ${LEAF_C}
check

# Link the Tcl extension as a shared library
echo "creating library..."
gcc -shared -o ${LEAF_SO} ${LEAF_O} ${USER_O} \
    -L ${CoMD_HOME}/src-flat -l CoMD_Shared \
    -Wl,-rpath -Wl,${CoMD_HOME}/src-flat
check
echo "created library: ${LEAF_SO}"

# Make the Tcl package index
export LEAF_PKG LEAF_SO LEAF_TCL LEAF_VERSION
${TCLSH} make-package.tcl > pkgIndex.tcl
check
echo "created package."

exit 0
