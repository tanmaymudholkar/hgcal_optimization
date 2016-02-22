#!/bin/bash

for f in *1000*; do
    echo $f
    newname=$(echo $f | sed 's|\(.*\)_sigreg1000\(.*\)|\1\2|')
    echo $newname
    mv $f $newname
done
