#!/bin/bash

cd userlib &&
    make clean &&
    make dictionary &&
    make -j 4 &&
    cd .. &&
    make &&
    echo "Build successful"
