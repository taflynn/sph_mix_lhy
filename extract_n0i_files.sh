#!/bin/bash

extracted_dir="./data/$1extracted"

echo $extracted_dir

mkdir $extracted_dir 

echo "script extracts:"
echo "1) del_n(r=0) [central density difference]"
echo "2) n_i(r=0) [centred by mean]"
echo "3) n_i(r=0) [uncentred]"
echo "4) n_i(r) [ground state wavefunction]"
echo "5) r [computation box in .npy format]"
echo "6) Theoretical parameters"
echo "7) Simulation input parameters"

for i in $(seq 1 1 $2);
do
  current_dir="$1$i"
  echo "Extracting data from directory:" $current_dir 
  mkdir "$extracted_dir/$current_dir"
  cp "./data/$current_dir/del_n0s.csv" "$extracted_dir/$current_dir/"
  cp "./data/$current_dir/n01_centred.csv" "$extracted_dir/$current_dir/"
  cp "./data/$current_dir/n02_centred.csv" "$extracted_dir/$current_dir/"
  cp "./data/$current_dir/n01_uncentred.csv" "$extracted_dir/$current_dir/"
  cp "./data/$current_dir/n02_uncentred.csv" "$extracted_dir/$current_dir/"
  cp "./data/$current_dir/imag_fin_dens1.csv" "$extracted_dir/$current_dir/"
  cp "./data/$current_dir/imag_fin_dens2.csv" "$extracted_dir/$current_dir/"
  cp "./data/$current_dir/r_array.npy" "$extracted_dir/$current_dir/"
  cp "./data/$current_dir/theory_params.json" "$extracted_dir/$current_dir/"
  cp "./data/$current_dir/config_dens_ulck$i.json" "$extracted_dir/$current_dir/"
done
