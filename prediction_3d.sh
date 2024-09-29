#! /bin/bash

k_clusters=(10000 20000 30000) #5000
m_clusters=(60 120 180 240 300 360 420)

for k in "${k_clusters[@]}"; do
    for m in "${m_clusters[@]}"; do
        ./bin/predict  --train_locs windspeed_3d/meta_train_1000000 --test_locs windspeed_3d/meta_test_1000000 --train_data windspeed_3d/observation_train_1000000 --test_data windspeed_3d/observation_test_1000000 -n 100000 -k $k -m $m --theta 0.43,0.00297189,0.45 --distance_metric 2 --num_threads 30 
    done
done

mkdir -p log/prediction_3d
mv log/conditional_simulation_k_* log/prediction_3d