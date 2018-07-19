#!/bin/bash

for rbind in .1 .09 .08 .07 .06 .05 .04 .03 .02 .01
do
	for koff in .1 .2 .5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000
	do
	tail -c 100000000 Data/Decay/rbind${rbind}/koff${koff}.txt > Data/Decay/rbind${rbind}/koff${koff}truncate.txt
	awk 'NR % 10 == 0' Data/Decay/rbind${rbind}/koff${koff}truncate.txt > Data/Decay/rbind${rbind}/koff${koff}trim.txt
	done 
done