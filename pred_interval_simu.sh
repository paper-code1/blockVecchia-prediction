#!/bin/bash

# iterate from 0 to 8, keep 6 digits
filename="params.txt"
k=200
m_all=(10 30 60 120 180 240 300)

# Loop through the lines in the file
while IFS=' ' read -r sigma beta nu seed; do
    echo "===================================================="
    echo "===================================================="
    echo "sigma: $sigma, beta: $beta, nu: $nu, seed:$seed"
    sigma_str=$(printf "%.6f" "$sigma")
    beta_str=$(printf "%.6f" "$beta")
    nu_str=$(printf "%.6f" "$nu")
    data_folder="./simu_ds/prediction_data_${beta_str}_${nu_str}"

    for i in {1..50}
    do
        for m in ${m_all[@]}
        do
            # file path using the variable i, e.g., LOC_1_train.csv, Z1_1_train.csv
            train_locs=$data_folder/LOC_${i}_train.csv
            test_locs=$data_folder/LOC_${i}_test.csv
            train_data=$data_folder/Z1_${i}_train.csv
            test_data=$data_folder/Z1_${i}_test.csv
            ./bin/predict  --train_locs $train_locs --test_locs $test_locs --train_data $train_data --test_data $test_data -n 2000 -k $k -m $m --theta ${sigma_str},${beta_str},${nu_str} --distance_metric 2 --num_threads 20 --seed $i
        done
    done

    mkdir -p ./log/pred_simu_${beta_str}_${nu_str}
    mv ./log/conditional_simulation* ./log/pred_simu_${beta_str}_${nu_str}
    echo "===================================================="
done < "$filename"