#!/bin/bash

name="aaa"
data="/data/anl2015/FEB2015/135SB/"

rm -f $name.* $name-tof*.dat

for i in `ls -tr $data/a135feb_12.ldf` 
#for i in `ls -tr $data/a135feb_12*.ldf` 
do 
#if [ "$i" == "$data/a135feb_12.ldf" ];
#then
#continue
#fi


     cmd=$cmd"file $i\ngo\n"

done
cmd=$cmd"end\n"
#echo -e $cmd
make -j16 && echo -e $cmd | ./pixie_ldf_c $name