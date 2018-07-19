#!/bin/bash
mkdir  /Users/allardlab/Smoldyn/Repots/SmoldynSimpleMFPT/Data/nSweep
for n in {2..18}
do
	./MFPT ${n} 3 0 .00001 1 100 Data/nSweep/${n}.txt
done