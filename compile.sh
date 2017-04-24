#!/bin/bash

nvcc -std=c++11 -Xptxas -g -G -arch=sm_20 -I /usr/local/cuda/include/ main_cuda.cpp input.cpp cudakernels.cu -o qs
