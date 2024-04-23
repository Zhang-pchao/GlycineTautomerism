#!/bin/bash

rm transitions.dat

touch transitions.dat

for (( i=101; i<161; i=i+1 )); do
	tail -n 1 $i/time >> transitions.dat	
done	
