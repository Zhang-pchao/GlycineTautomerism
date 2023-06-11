#!/bin/bash

for x in run{0..20}
do
tail -1 $x/time | awk '{print $5}' >> times.dat
done
