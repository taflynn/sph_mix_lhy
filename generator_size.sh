#!/bin/bash
x=1
let "N_steps = ($3 - $1) / $2"
sig_step=$(bc<<<"scale=2; ($5 - $4) / $N_steps")
let "gauss_sig = $4"
for i in $(seq $1 $2 $3);
do
    cp ./config_dens_ulck.json ./config_dens_ulck$x.json
    touch ./config_dens_ulcktemp1.json
    touch ./config_dens_ulcktemp2.json
    awk -v val=$gauss_sig '{ gsub("\"GAUSS_SIGMA\":[0-9]+.0","\"GAUSS_SIGMA\":"val,$0);print $0}' ./config_dens_ulck.json > ./config_dens_ulcktemp1.json
    let "imb = $i / $6 + $i"
    awk -v val=$imb '{ gsub("\"N1\":[0-9]+","\"N1\":"val,$0);print $0}' ./config_dens_ulcktemp1.json > ./config_dens_ulcktemp2.json
    awk -v val=$i '{ gsub("\"N2\":[0-9]+","\"N2\":"val,$0);print $0}' ./config_dens_ulcktemp2.json > ./config_dens_ulck$x.json
    ((x = x + 1)) 
    gauss_sig=$(bc <<<"scale=2; $gauss_sig + $sig_step")
done
rm ./config_dens_ulcktemp1.json
rm ./config_dens_ulcktemp2.json
