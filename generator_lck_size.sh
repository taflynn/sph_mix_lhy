#!/bin/bash
x=1
let "N_steps = ($3 - $1) / $2"
sig_step=$(bc<<<"scale=2; ($5 - $4) / $N_steps")
let "gauss_sig = $4"
touch ./config_dens_lcktemp.json
for i in $(seq $1 $2 $3);
do
    cp ./config_dens_lck.json ./config_dens_lck$x.json
    awk -v val=$gauss_sig '{ gsub("\"GAUSS_SIGMA\":[0-9]+.0","\"GAUSS_SIGMA\":"val,$0);print $0}' ./config_dens_lck.json > ./config_dens_lcktemp.json
    awk -v val=$i '{ gsub("\"N\":[0-9]+","\"N\":"val,$0);print $0}' ./config_dens_lcktemp.json > ./config_dens_lck$x.json
    ((x = x + 1)) 
    gauss_sig=$(bc <<<"scale=2; $gauss_sig + $sig_step")
done
rm config_dens_lcktemp.json
