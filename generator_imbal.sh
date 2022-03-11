#!/bin/bash
x=1
for i in $1 $2 $3 $4 $5
do
    cp ./config_dens_ulck.json ./config_dens_ulck$x.json
    awk -v val=$i '{ gsub("\"N1\":[0-9]+","\"N1\":"val,$0);print $0}' ./config_dens_ulck.json > ./config_dens_ulck$x.json
    ((x = x + 1)) 
done
