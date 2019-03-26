#!/bin/bash
set -e

here=`pwd`

case "`uname`" in
  Darwin*)
    sedi='sed -i ""'
    setsedi='s|sed -i -e|sed -i "" -e|'
    ;;
  *)
    sedi='sed -i'
    setsedi='s|sed -i "" -e|sed -i -e|'
    ;;
esac

$sedi -e"$setsedi" create.sh
$sedi -e"$setsedi" example/run.sh 
$sedi -e"$setsedi" example_cpp/run.sh 
$sedi -e"$setsedi" checks/run.sh

sedArg="s| *amp4hefDir *=.*|amp4hefDir=$here|"
$sedi -e"$sedArg" create.sh
$sedi -e"$sedArg" example/run.sh
$sedi -e"$sedArg" example_cpp/run.sh
$sedi -e"$sedArg" checks/run.sh
$sedi -e"$sedArg" applications/two2two/compile.sh

