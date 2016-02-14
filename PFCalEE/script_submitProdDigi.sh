#!/bin/bash

for i in {1..100}; do
    sleep 2
    ./submitProdDigi.py -s 2nd -q 1nw -t V05-02-04 -v 34 -m 2 -a 2.1 -b 0 -d gamma -n 1000 -g -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/vcnncal -R $i
done


# ./submitProdDigi.py -s 2nd -q 1nw -t V05-02-04 -v 40 -m 2 -a 2.1 -b 0 -d gamma -n 1000 -g -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/vfrac1
# ./submitProdDigi.py -s 2nd -q 1nw -t V05-02-04 -v 41 -m 2 -a 2.1 -b 0 -d gamma -n 1000 -g -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/vfrac75
# ./submitProdDigi.py -s 2nd -q 1nw -t V05-02-04 -v 42 -m 2 -a 2.1 -b 0 -d gamma -n 1000 -g -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/vfrac50
# ./submitProdDigi.py -s 2nd -q 1nw -t V05-02-04 -v 43 -m 2 -a 2.1 -b 0 -d gamma -n 1000 -g -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/vfrac25
# ./submitProdDigi.py -s 2nd -q 1nw -t V05-02-04 -v 44 -m 2 -a 2.1 -b 0 -d gamma -n 1000 -g -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/vflat
