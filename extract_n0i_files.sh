#!/bin/bash

extracted_dir="./data/$1extracted"

echo $extracted_dir

mkdir $extracted_dir 

for i in $(seq 1 1 $2);
do
  current_dir="$1$i"
  echo $current_dir 
  mkdir "$extracted_dir/$current_dir"
  cp "./data/$current_dir/n01.csv" "$extracted_dir/$current_dir/"
  cp "./data/$current_dir/n02.csv" "$extracted_dir/$current_dir/"
  cp "./data/$current_dir/imag_fin_dens1.csv" "$extracted_dir/$current_dir/"
  cp "./data/$current_dir/imag_fin_dens2.csv" "$extracted_dir/$current_dir/"
  cp "./data/$current_dir/theory_params.json" "$extracted_dir/$current_dir/"
  cp "./data/$current_dir/config_dens_ulck$i.json" "$extracted_dir/$current_dir/"
done
