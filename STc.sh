#!/bin/bash -e
a=1
c=1
while [ $a -lt 12 ]
do 
  echo $a
  b=1
  while [ $b -lt 12 ]
  do 
    echo $b
    gfortran -cpp -Dpyweath -Dvar1=${a} -Dvar2=${b} -o run${c} o2profile+silweath+o2_v9_7.f -lopenblas -g -fcheck=all 
    cp run${c} ./run
    rm run${c}
    b=`expr $b + 1`
    c=`expr $c + 1`
  done
  a=`expr $a + 1`
done
