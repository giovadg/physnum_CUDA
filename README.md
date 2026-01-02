#  physnum_CUDA
C++ GPU/CUDA accelerated numerical algorithms.

This repository contains some of the scripts developed for the Physique-Numerique class held by Prof. Laurent Villard at EPFL.

The numerical algorithms are here GPU-accelerated (CPU version only has been developed in the framework of the Physique numerique class) with the fewest modifications compared to the CPU version.

## Requirements

- The MonteCarlo algorithm in this project requires the boost library for random number generation 
- CUDA libreries
  
## Benchmark
To test the GPU acceleration there is the comparison with the CPU baseline. Possible by using the correct input parameter and run the algorithms on host -still relying on CUDA libreries for comparison but without the acceleration -

## Format
The codes are Jupp. nootebok (.ipynb) that can be open and run on colab to exploit the acceleration.

## Output
The notebooks compile and run. Finally some physical quantities are plot for comparison.

