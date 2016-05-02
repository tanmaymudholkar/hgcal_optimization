#!/usr/bin/env python

# import os,sys
# import optparse
# import commands
# import math
# import random

# random.seed()


# version 34

twcu_t = 39.6
tw_t = 41.6
x0wcu = 5.1223
x0w = 3.50418
rwcu_t = twcu_t/x0wcu
rw_t = tw_t/x0w

x0neutmod = 471.303
x0al = 88.9632
x0foam = 4482.29
x0air = 303921
x0pcb = 187.31
x0si = 93.6607
x0cu = 14.3558
x0cfm = 3494.61
lneutmod = 100
lal = 4.0
lfoam = 26.0
lcueven = 1.0
lcfm = 1.0
lair = 3.0
lpcb = 2.0
lsi = 0.3
lcuodd = 6.0;

rw_min = lcueven/x0cu + lcfm/x0cfm + lair/x0air + lpcb/x0pcb + lsi/x0si
rwcu_min = lcuodd/x0cu + lsi/x0si + lpcb/x0pcb + lair/x0air

totrlen = 0
weights = []
# layer 0
totrlen = lneutmod/x0neutmod + lal/x0al + lfoam/x0foam + lcfm/x0cfm + 2.6/x0w + 0.5/x0cu + lair/x0air + lpcb/x0pcb + lsi/x0si
weights.append(totrlen)
# print "cumulative rlen ",totrlen

def rw(lw):
    return rw_min + lw/x0w

def rwcu(lwcu):
    return rwcu_min + 2*lwcu/x0wcu

# layer 1
totrlen += rwcu(1)
weights.append(rwcu(1))
# print "cumulative rlen ",totrlen

# layers 2 ---> 7
for i in range(0,3):
    totrlen += rw(2.6)
    weights.append(rw(2.6))
    # print "cumulative rlen ",totrlen
    totrlen += rwcu(1)
    weights.append(rwcu(1))
    # print "cumulative rlen ",totrlen

# layers 8 ---> 15
for i in range(0,4):
    totrlen += rw(3.6)
    weights.append(rw(3.6))
    # print "cumulative rlen ",totrlen
    totrlen += rwcu(1.75)
    weights.append(rwcu(1.75))
    # print "cumulative rlen ",totrlen

# layers 16 ---> 23
for i in range(0,4):
    totrlen += rw(4.2)
    weights.append(rw(4.2))
    # print "cumulative rlen ",totrlen
    totrlen += rwcu(2.2)
    weights.append(rwcu(2.2))
    # print "cumulative rlen ",totrlen


print "version 34:"
print "_________________________________________________________________________________"
print weights
normalized_weights=[1.0]
for i in range(1,len(weights)):
    normalized_weights.append(weights[i]/weights[0])

for i in range(0,len(normalized_weights)):
    print "%.3f"%(normalized_weights[i])
    
print "cumulative rlen ",totrlen
print "_________________________________________________________________________________"


# version 30

twcu_t = 39.6
tw_t = 41.6
x0wcu = 5.1223
x0w = 3.50418
rwcu_t = twcu_t/x0wcu
rw_t = tw_t/x0w

x0neutmod = 471.303
x0al = 88.9632
x0foam = 4482.29
x0air = 303921
x0pcb = 187.31
x0si = 93.6607
x0cu = 14.3558
x0cfm = 3494.61
lneutmod = 100
lal = 4.0
lfoam = 26.0
lcueven = 1.0
lcfm = 1.0
lair = 3.0
lpcb = 2.0
lsi = 0.3
lcuodd = 6.0;

rw_min = lcueven/x0cu + lcfm/x0cfm + lair/x0air + lpcb/x0pcb + lsi/x0si
rwcu_min = lcuodd/x0cu + lsi/x0si + lpcb/x0pcb + lair/x0air

totrlen = 0
weights = []
# layer 0
totrlen = lneutmod/x0neutmod + lal/x0al + lfoam/x0foam + lcfm/x0cfm + 2.0/x0w + 0.5/x0cu + lair/x0air + lpcb/x0pcb + lsi/x0si
weights.append(totrlen)
# print "cumulative rlen ",totrlen

def rw(lw):
    return rw_min + lw/x0w

def rwcu(lwcu):
    return rwcu_min + 2*lwcu/x0wcu

# layer 1
totrlen += rwcu(0.6)
weights.append(rwcu(0.6))
# print "cumulative rlen ",totrlen

# layers 2 ---> 9
for i in range(0,4):
    totrlen += rw(2.0)
    weights.append(rw(2.0))
    # print "cumulative rlen ",totrlen
    totrlen += rwcu(0.6)
    weights.append(rwcu(0.6))
    # print "cumulative rlen ",totrlen

# layers 10 ---> 19
for i in range(0,5):
    totrlen += rw(2.8)
    weights.append(rw(2.8))
    # print "cumulative rlen ",totrlen
    totrlen += rwcu(1.2)
    weights.append(rwcu(1.2))
    # print "cumulative rlen ",totrlen

# layers 20 ---> 27
for i in range(0,4):
    totrlen += rw(4.2)
    weights.append(rw(4.2))
    # print "cumulative rlen ",totrlen
    totrlen += rwcu(2.2)
    weights.append(rwcu(2.2))
    # print "cumulative rlen ",totrlen


print "version 30:"
print "_________________________________________________________________________________"
print weights
normalized_weights=[1.0]
for i in range(1,len(weights)):
    normalized_weights.append(weights[i]/weights[0])

for i in range(0,len(normalized_weights)):
    print "%.3f"%(normalized_weights[i])
    
print "cumulative rlen ",totrlen
print "_________________________________________________________________________________"
