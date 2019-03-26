#!/bin/bash
set -e

amp4hefDir=/home/user0/repos/amp4hef
srcDir=$amp4hefDir/src
buildDir=$amp4hefDir/build
compiler="gfortran"

mkdir -p $buildDir

function info {
  echo ""
  echo "======================================================================="
  echo " Execute"
  echo ""
  echo " $ $0 source"
  echo ""
  echo " to create the source file"
  echo ""
  echo " $buildDir/amp4hef.f03"
  echo ""
  echo " which you can compile to make the Fortran  module amp4hef  available."
  echo " To create a library  libamp4hef.a  in the build directory, execute"
  echo ""
  echo " $ $0 lib"
  echo ""
  echo " The gfortran compiler is used in this case. Alternatively, execute"
  echo ""
  echo " $ $0 lib -compiler 'yourFortranCompiler -yourFlags'"
  echo ""
  echo " To choose another build directory, add the flag"
  echo ""
  echo " -build yourBuildDirectory"
  echo ""
  echo " It is preferable to use this flag instead of copying files, because"
  echo " the necessary paths are set automatically."
  echo " The tasks source, lib etc. can also be indicated via a flag, eg."
  echo ""
  echo " $ $0 \\"
  echo "     -task lib \\"
  echo "     -compiler 'yourFortranCompiler -yourFlags' \\"
  echo "     -build yourBuildDirectory"
  echo ""
  echo "======================================================================="
  echo ""
}

if [ "$#" -eq 0 ]; then
  info
  exit 0
fi

function create_source {
  cat \
    $srcDir/amp4hef_io.f03 \
    $srcDir/amp4hef_qomentum.f03 \
    $srcDir/amp4hef_aux.f03 \
    $srcDir/amp4hef_ng.f03 \
    $srcDir/amp4hef_qq.f03 \
    $srcDir/amp4hef_main.f03 \
    > $buildDir/amp4hef.f03
  sed -i -e"s|\!(path_tbldir\!)|'$buildDir/'|" $buildDir/amp4hef.f03
  cp $srcDir/amp4hef.tbl $buildDir
}

args=("$@")
i="0"
while [ $i -lt $# ]; do
  case "${args[$i]}" in
  "-h"|"-help"|"--help"|"-i"|"-info"|"--info")
    info
    exit 0
    ;;
  "-task")
    i=$[$i+1]
    task=${args[$i]}
    ;;
  "-build")
    i=$[$i+1]
    buildDir=${args[$i]}
    ;;
  "-compiler")
    i=$[$i+1]
    compiler=${args[$i]}
    ;;
  *)
    if [ $i -eq 0 ]; then
      task=${args[$i]}
    else
      echo "ERROR: option ${args[$i]} not defined"
      exit 1
    fi
    ;;
  esac
  i=$[$i+1]
done

case $task in

"source")
  create_source
  ;;

"lib")
  create_source
  cd $buildDir
  $compiler -c amp4hef.f03
  ar cru libamp4hef.a amp4hef.o
  ranlib libamp4hef.a
  ;;

"clean")
  rm -r -f $buildDir/*
  rm -r -f $amp4hefDir/example/build/*
  rm -r -f $amp4hefDir/example_cpp/build/*
  rm -r -f $amp4hefDir/checks/build/*
  ;;

*)
  info
  ;;

esac

