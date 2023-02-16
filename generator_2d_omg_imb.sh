#!/bin/bash
x=1
for j in $(seq $1 $2 $3);
do
  for i in $(seq $4 $5 $6);
  do
    cp ./config_dens_ulck.json ./config_dens_ulck$x.json
    touch ./config_dens_ulcktemp1.json
    touch ./config_dens_ulcktemp2.json
    awk -v val=$j '{ gsub("\"OMEGA1\":[0-9]+.0","\"OMEGA1\":"val,$0);print $0}' ./config_dens_ulck.json > ./config_dens_ulcktemp1.json
    awk -v val=$j '{ gsub("\"OMEGA2\":[0-9]+.0","\"OMEGA2\":"val,$0);print $0}' ./config_dens_ulcktemp1.json > ./config_dens_ulcktemp2.json
    awk -v val=$i '{ gsub("\"N1\":[0-9]+","\"N1\":"val,$0);print $0}' ./config_dens_ulcktemp2.json > ./config_dens_ulck$x.json
    ((x = x + 1)) 
  done
done
rm ./config_dens_ulcktemp1.json
rm ./config_dens_ulcktemp2.json
