#!/bin/bash
x=1
for i in $(seq $1 $2 $3);
do
    cp ./config_dens_ulck.json ./config_dens_ulck$x.json
    awk -v val=$i '{ gsub("\"OMEGA2\":[0-9]+","\"OMEGA2\":"val,$0);print $0}' ./config_dens_ulck.json > ./config_dens_ulck$x.json
    ((x = x + 1)) 
done
