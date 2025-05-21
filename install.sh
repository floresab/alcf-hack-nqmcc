#!/bin/bash

set -e # exit on error

usage()
{
  echo "nQMCC Install Script"
  echo ""
  echo "USAGE:"
  echo "./install-nqmcc.sh [--help|-h] | [-p|--prefix=PREFIX] | [-m|--machine=MACHINE] | [-g|--debug=DEBUG] | [-v|--VERBOSE=VERBOSE] "
  echo ""
  echo " --help             Display this help text"
  echo " --prefix=PREFIX    Install binary in 'PREFIX/bin'"
  echo "                    Default prefix='Path/To/nQMCC/build/bin'"
  echo " --debug=DEBUG      Compile with -g"
  echo "                    Default off"
  echo " --verbose=VERBOSE  Print Compiler information"
  echo "                    Default off"
  echo " --machine=MACHINE  sets compiler and flags"
  echo "                    Default MACHINE='gnu'"
  echo "                    Available:          '"
  echo "                               gnu-cpu  '"
  echo "                               gnu-gpu'"
  echo "                               improv'"
  echo "                               aurora'"
  echo "                               nersc-cpu'"
  echo "                               nersc-gpu'"
#  echo "                               polaris'"
  echo ""
}

PREFIX="$(pwd)/build"
MACHINE="gnu-cpu"
DEBUG=""
VERBOSE=""
PHI_R4=0
while [ "$1" != "" ]; do
  PARAM=$(echo "$1" | awk -F= '{print $1}')
  VALUE=$(echo "$1" | awk -F= '{print $2}')
  case $PARAM in
    -h | --help)
      usage
      exit
      ;;
    -p | --prefix)
      PREFIX=$VALUE
      ;;
    -m | --machine)
      MACHINE=$VALUE
      ;;
    -g | --debug)
      DEBUG=" -g "
      ;;
    -v | --verbose)
      VERBOSE=" --verbose"
      ;;
    *)
      echo "ERROR: unknown parameter \"$PARAM\""
      usage
      exit 1
      ;;
  esac
  shift
done

set -u # error on use of undefined variable

case $MACHINE in
  gnu-cpu)
  CMK_FC=${CMK_FC:-"mpif90"}
  #CMK_FC=${CMK_FC:-"gfortran-14"}
  CMK_FLAG="-Wall -Dgpu_offload=0 -cpp -O3 -fopenmp -std=f2018 -ffree-line-length-512 -fbacktrace -fcheck=all -pedantic"
  ;;
  aurora)
  CMK_FC=${CMK_FC:-"mpif90"}
  CMK_FLAG="-warn -Dgpu_offload=1 -cpp -O3 -fiopenmp -fopenmp-targets=spir64_gen -Xopenmp-target-backend=spir64_gen \"-device pvc\" "
  ;;
  nersc)
  CMK_FC=${CMK_FC:-"ftn"}
  CMK_FLAG="-Wall -Dgpu_offload=1 -fast -O3 -mp=gpu,multicore -cpp -gpu=cc80 -Minfo=mp,accel"
  ;;

  *)
  echo "ERROR: unknown machine \"$MACHINE\""
      usage
      exit 1
      ;;
esac

CI=${CI:-"false"} # GitHub Actions workflows set CI=true
CMK_FLAG=${CMK_FLAG}" ${DEBUG} ${VERBOSE}"

#buildit
cmake -S . -B ${PREFIX} -DCMAKE_Fortran_COMPILER=${CMK_FC}  \
                        -DCMAKE_Fortran_FLAGS="${CMK_FLAG}"\
