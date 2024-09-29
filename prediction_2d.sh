#! /bin/bash

k_clusters=(5000 10000 20000 30000)
m_clusters=(60 120 180 240 300 360 420)

for k in "${k_clusters[@]}"; do
    for m in "${m_clusters[@]}"; do
        ./bin/predict  --train_locs ./soil-moi/meta_train_0.9 --test_locs ./soil-moi/meta_test_0.9 --train_data ./soil-moi/observation_train_0.9 --test_data ./soil-moi/observation_test_0.9 -n 200000 -k $k -m $m --scale_factor 9348.317 --theta 0.65,0.00128365,0.475 --distance_metric 1 --num_threads 15
    done
done

mkdir -p log/prediction_2d
mv log/conditional_simulation_k_* log/prediction_2d