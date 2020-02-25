#!/bin/bash

mkdir distFeats
for f in *.tsv *.txt; do
    #open the file and get the number of trees
    echo $f
    line=$(sed -n 2p "$f")
    echo $line
    NUMBER=$(echo $line | grep -o -E '[0-9]+')
  

    let NUMBER-=1
    #loop through for each tree and find the the distinguishing features
    until [  $NUMBER -lt 0 ]; do
             echo COUNTER $NUMBER
             generateDFF $f -i $NUMBER
             let NUMBER-=1
    done
done
