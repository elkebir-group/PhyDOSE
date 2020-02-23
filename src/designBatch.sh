#!/bin/bash

mkdir distFeats
for f in ./*; do
    #open the file and get the number of trees
    echo $f
    line=$(sed -n '2{p;q}' $f)
    echo $line
    NUMBER=$(echo $line | grep -o -E '[0-9]+')
  

    
    #loop through for each tree and find the the distinguishing features
    until [  $NUMBER -lt 1 ]; do
             echo COUNTER $NUMBER
             design $f -i $NUMBER
             let NUMBER-=1
    done
done
