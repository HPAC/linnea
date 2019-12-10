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
echo "cmake v3.13.2"
echo "clang 8.0"
echo "python 3.6"
echo "Intel MKL 2019"

echo "Linnea environment dir: $LINNEA_INSTALL_DIR"
echo
echo "############################################"
echo

module purge
module load DEVELOP
module load LIBRARIES
module load clang/8.0 && CLANG_MODULE="YES" || CLANG_MODULE="NO"
module load python/3.6.0 && PYTHON_MODULE="YES" || PYTHON_MODULE="NO"
module load cmake/3.13.2 && CMAKE_MODULE="YES" || CMAKE_MODULE="NO"
module load intelmkl/2019 && INTEL_MKL_MODULE="YES" || INTEL_MKL_MODULE="NO"

echo "############################################"
echo "Python 3.6.0:..." $PYTHON_MODULE
echo "CMake 3.13.2:..." $CMAKE_MODULE
echo "Intel MKL 2019:....." $INTEL_MKL_MODULE
echo "Clang 8.0:.........." $CLANG_MODULE
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

echo "Installing Eigen"
    cd $SRC_DIR/eigen
    git checkout 3.3.7
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$LIB_DIR" ..
    make install

echo "Installing Armadillo"
    cd $SRC_DIR/armadillo
    git checkout 9.800.x
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
    python3.6 -m venv $VENV

echo "Installing Linnea in the virtual environment"
    source $VENV/bin/activate
    cd $SRC_DIR/linnea
    pip install --upgrade pip setuptools
    pip install -e .
    deactivate

echo "Installing Julia"
    cd $SRC_DIR/julia
    git checkout v1.3.0
    echo "USE_INTEL_MKL = 1" > $SRC_DIR/julia/Make.user
    make -j 24
    ./julia -e "using Pkg; Pkg.add(PackageSpec(url=\"https://github.com/HPAC/MatrixGenerator.jl.git\", rev=\"master\"))"

echo "Finished."

