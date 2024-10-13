# Block Vecchia Prediction

## 1.Introduction

Block Vecchia Prediction is a powerful tool for spatial prediction and conditional simulation using the Block Vecchia approximation method. This package is designed for efficient handling of large-scale spatial datasets, particularly useful in geostatistics, environmental science, and other fields dealing with spatially correlated data.

The implementation uses C++ for high-performance computing and supports multi-threading for improved efficiency on modern multi-core processors.

## 2.Installation

### 2.1 Prerequisites

- C++ compiler with C++17 support (e.g., g++ 10.2.0 or later)
- OpenMP (for multi-threading support 4.1.0 or later)
- GSL (GNU Scientific Library 2.6.0 or later)
- BLAS (Basic Linear Algebra Subprograms)

### 2.2 Building from source


1. Ensure you have the required libraries installed. On Ubuntu or Debian-based systems, you can install them with:
   ```
   sudo apt-get install libgsl-dev libblas-dev libopenmpi-dev
   ```

This will create the `predict` executable in the `bin` directory.

## 3.Usage

The main executable `predict` can be run with various command-line arguments to specify input data, parameters, and output options.

### 3.1 Basic usage:

```
bash
./bin/predict --train_locs <train_locations_file> --test_locs <test_locations_file> \
--train_data <train_data_file> --test_data <test_data_file> \
-n <num_samples> -k <num_clusters> -m <block_size> \
--scale_factor <scale> --theta <theta_params> \
--distance_metric <metric> --num_threads <threads> --seed <seed>
```

This script will run the prediction for various combinations of k (number of clusters) and m (block size) values, and move the output logs to the `log/prediction_2d` directory.

## Output

The program will generate output files in the `log` directory, with filenames indicating the parameters used for each run.

## License

Apache License 2.0