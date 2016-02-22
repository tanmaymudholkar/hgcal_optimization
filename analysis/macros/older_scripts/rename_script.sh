#!/bin/bash

for file in plot_*; do
    newname=$(echo $file | sed 's|plot_v40|plot_compare_v40|')
    mv $file $newname
done
