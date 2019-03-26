#!/bin/bash
set -e

FC="gfortran"
#FC="nagfor"
#FC="ifort"

here=`pwd`
amp4hefDir=/home/user0/repos/amp4hef
buildDir=$here/build
mkdir -p $buildDir

$amp4hefDir/create.sh -task source -build $buildDir
cp $here/src/main.f03 $buildDir

cd $buildDir
$FC amp4hef.f03 main.f03
mv a.out $here/main.out
cd $here

