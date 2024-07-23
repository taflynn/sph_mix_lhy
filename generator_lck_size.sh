#!/bin/bash
x=1
for i in $(seq $1 $2 $3);
do
    cp ./config_dens_lck.json ./config_dens_lck$x.json
    awk -v val=$i '{ gsub("\"N\":[0-9]+","\"N\":"val,$0);print $0}' ./config_dens_lck.json > ./config_dens_lck$x.json
    ((x = x + 1)) 
done
