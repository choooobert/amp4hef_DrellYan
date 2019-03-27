#!/bin/bash
set -e

FC="gfortran"
here=`pwd`
amp4hefDir=/home/hubert/repo/amp4hef_DrellYan
exampleDir=$amp4hefDir/example_cpp
buildDir=$exampleDir/build
mkdir -p $buildDir

if [ "$#" -eq 0 ]; then
   echo ""
   echo "================================================================"
   echo " Execution:"
   echo ""
   echo " $ $0 ../example/datafiles/gs_gs_3g.evt"
   echo ""
   echo " Instead of gs_gs_3g.evt, you can choose any file in the"
   echo " directory datafiles."
   echo ""
   echo " Cleaning:"
   echo ""
   echo "$ $0 clean"
   echo "================================================================"
   echo ""
   exit 0
fi

if [  "$1" = "clean" ]; then
   rm -r -f $buildDir/*
   $amp4hefDir/create.sh clean
   exit 0
fi

$amp4hefDir/create.sh -task source -build $buildDir
cp $exampleDir/src/main.cpp $buildDir
cp $exampleDir/src/amp4hef.hpp $buildDir

cd $buildDir
gfortran -c amp4hef.f03 -o amp4hef.o
g++  -g -std=c++11 -Wall -c main.cpp -o main.o
g++ -Wall main.o amp4hef.o  -o program.out -lm -lgfortran
cd $here

time $buildDir/program.out $1

