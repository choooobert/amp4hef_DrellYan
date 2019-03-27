#!/bin/bash
set -e

FC="gfortran"
#FC="nagfor"
#FC="ifort"

here=`pwd`
amp4hefDir=/home/hubert/repo/amp4hef_DrellYan
exampleDir=$amp4hefDir/example
buildDir=$exampleDir/build
mkdir -p $buildDir

if [ "$#" -eq 0 ]; then
   echo ""
   echo "================================================================"
   echo " Choose your Fortran compiler at the beginning of the file"
   echo ""
   echo " $0"
   echo ""
   echo " and execute"
   echo ""
   echo " $ $0 datafiles/gs_gs_3g.evt"
   echo ""
   echo " Instead of gs_gs_3g.evt, you can choose any file in the"
   echo " directory datafiles."
   echo "================================================================"
   echo ""
   exit 0
fi

$amp4hefDir/create.sh -task source -build $buildDir
cp $exampleDir/src/main.f03 $buildDir

cd $buildDir
$FC amp4hef.f03 main.f03
cd $here

time $buildDir/a.out $1

