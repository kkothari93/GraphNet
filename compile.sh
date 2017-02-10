#!/bin/bash

nvcc -std=c++11 -arch=sm_30 -I /usr/local/cuda/include/ main_cuda.cpp input.cpp cudakernels.cu -o qs
