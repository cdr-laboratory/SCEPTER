#!/bin/bash -e
a=1
c=1
while [ $a -lt 12 ]
do 
  # echo $a
  b=1
  while [ $b -lt 12 ]
  do 
    # echo $b
    echo $c
    # if [ $c -eq 10 ]||[ $c -eq 15 ]||[ $c -eq 16 ]||[ $c -eq 14 ]||[ $c -eq 18 ]||[ $c -eq 19 ]||[ $c -eq 27 ]; then 
      # b=`expr $b + 1`
      # c=`expr $c + 1`
      # continue
    # fi
    gfortran -cpp -Dpyweath -Dvar1=${a} -Dvar2=${b} -o run${c} o2profile+silweath+o2_v9_7.f -lopenblas -g -fcheck=all 
    if [[ -f ~/./run/run${c} ]]; then
      rm ./run/run${c}
    fi
    cp run${c} ./run
    rm run${c}
    b=`expr $b + 1`
    c=`expr $c + 1`
  done
  a=`expr $a + 1`
done
