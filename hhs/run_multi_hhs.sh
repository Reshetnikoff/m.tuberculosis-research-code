#!/bin/bash
# Script for running HHS programm for each drug and for each fold

# Configutation part #

# Directory with HHS format files
input_dir=../db/cv_bess_files
# Output directory (create three dir: {output_dir}/1, {output_dir}/3, {output_dir}/5
# before start)
output_dir=./output_hhs
# Hyperparameter - minimum frequency of mutations that lead to the resistance phenotype
filter=5

######################

drugs=(Kanamycin Capreomycin Ciprofloxacin Moxifloxacin Amikacin Prothionamide Ethionamide Ofloxacin Para-aminosalisylic_acid Isoniazid Ethambutol Pyrazinamide Rifampicin Streptomycin)
inds=(0 1 2 3 4)

for drug in ${drugs[@]}
do
for i in ${inds[@]}
do
../../programms/HHS/target/release/hhs -g "${input_dir}/${drug}.gt.domains.3.train.${i}.gz" -p "${input_dir}/${drug}.phen.domains.3.train.${i}" -d "${input_dir}/${drug}.gt.domains.3.train.${i}.dist.gz"  -t 6 --p1g1_filter $filter -o "${output_dir}/$filter/${drug}.gt.domains.3.train.${i}.result"
sleep 1
done
done

