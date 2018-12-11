#!/bin/bash

LINNEA_INSTALL_DIR=$HOME/Linnea
LIB_DIR=$LINNEA_INSTALL_DIR/lib
SRC_DIR=$LINNEA_INSTALL_DIR/src
VENV=$LINNEA_INSTALL_DIR/linnea.venv

PROMPT=1
set -e

function prompt {
  if [ "$PROMPT" -eq 1 ]; then
    answer=n
    while [ "$answer" != "Y" ]; do
      echo -n "Proceed with installation? (Y/n): "
      read answer
      if [ "$answer" = "n" ]; then
        echo "Quiting Installation..."
        exit 0
      fi
    done
  fi
}

echo "############################################"
echo "Please make sure you have the following packages installed:"
echo "cmake v3.10.1"
echo "gcc v7"
echo "python 3.6"
echo "Intel MKL"

echo "Linnea environment dir: $LINNEA_INSTALL_DIR"
echo
echo "############################################"
echo

module load python/3.6.0 && PYTHON_MODULE="YES" || PYTHON_MODULE="NO"
module load cmake/3.10.1 && CMAKE_MODULE="YES" || CMAKE_MODULE="NO"
module switch intel intel/18.0 && INTEL_MODULE="YES" || INTEL_MODULE="NO"
module load gcc/7 && GCC_MODULE="YES" || GCC_MODULE="NO"

echo "############################################"
echo "Python 3.6.0:..." $PYTHON_MODULE
echo "CMake 3.10.1:..." $CMAKE_MODULE
echo "Intel 18.0:....." $INTEL_MODULE
echo "Gcc 7:.........." $GCC_MODULE
echo "############################################"
echo "If all above modules are loaded you should proceed with the installation."
prompt

if [ -d $LINNEA_INSTALL_DIR ]; then
  rm -rf $LINNEA_INSTALL_DIR
fi

mkdir -p $LIB_DIR
mkdir -p $SRC_DIR
mkdir -p $VENV
echo
echo "Cloning all sources to $SRC_DIR ..."
echo
git clone https://github.com/HPAC/MatrixGenerator.jl.git $SRC_DIR/Matrixgenerator.jl
git clone https://github.com/HPAC/MatrixGeneratorMatlab.git $LIB_DIR/MatrixGeneratorMatlab
git clone https://github.com/HPAC/MatrixGeneratorCpp.git $SRC_DIR/MatrixGeneratorCpp
git clone https://github.com/eigenteam/eigen-git-mirror.git $SRC_DIR/eigen
git clone https://gitlab.com/conradsnicta/armadillo-code.git $SRC_DIR/armadillo
git clone https://github.com/JuliaLang/julia.git $SRC_DIR/julia
git clone https://github.com/HPAC/linnea.git $SRC_DIR/linnea

echo "Installing Julia" 
    cd $SRC_DIR/julia 
    git checkout 1c6f89f04a1ee4eba8380419a2b01426e84f52aa # Julia 1.1.0-DEV.468 from October 17, 2018
    echo "USE_INTEL_MKL = 1" > $SRC_DIR/julia/Make.user 
    make -j 1
    ./julia -e "using Pkg; Pkg.add(PackageSpec(url=\"https://github.com/HPAC/MatrixGenerator.jl.git\", rev=\"master\"))" 

echo "Installing Eigen" 
    cd $SRC_DIR/eigen 
    git checkout 3.3.5 
    mkdir build 
    cd build 
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$LIB_DIR" .. 
    make install 

echo "Installing Armadillo" 
    cd $SRC_DIR/armadillo 
    git checkout 9.200.x 
    mkdir build 
    cd build 
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$LIB_DIR" .. 
    make install 

echo "Installing MatrixGeneratorCpp" 
    cd $SRC_DIR/MatrixGeneratorCpp 
    mkdir build 
    cd build 
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$LIB_DIR" -DCMAKE_PREFIX_PATH="$LIB_DIR" -DWITH_ARMADILLO=On -DWITH_EIGEN=On -DWITH_MKL=On .. 
    make install 

echo "Creating Linnea virtual environment" 
    python3.6 -m pip install virtualenv --user 
    python3.6 -m virtualenv $VENV --always-copy 

echo "Installing Linnea in the virtual environment" 
    source $VENV/bin/activate 
    cd $SRC_DIR/linnea 
    pip install -e . 
    deactivate 

echo "Finished."

