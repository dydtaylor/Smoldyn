#!/bin/bash

for rbind in 1 .5 .4 .3 .2 .1 .09 .08 .07 .06 .05 .04 .03 .02 .01 .005 .002 .001
do
	mkdir  /Users/allardlab/Smoldyn/Repots/SmoldynDimerization/Data/Equilibrium/rbind${rbind};
done


for rbind in 1 .5 .4 .3 .2 .1 .09 .08 .07 .06 .05 .04 .03 .02 .01 .005 .002 .001
do
	for koff in .1 .2 .5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000
	do
	./dimer .00001 12.5 10 1000 ${koff} ${rbind} 0 Data/Equilibrium/rbind${rbind}/koff${koff}.txt
	done
done