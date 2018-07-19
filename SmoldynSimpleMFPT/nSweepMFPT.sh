#!/bin/bash
mkdir  /Users/allardlab/Smoldyn/Repots/SmoldynSimpleMFPT/Data/nSweep
for n in {2..18}
do
	./MFPT ${n} 100 3 0 Data/nSweep/${n}.txt
done