#!/bin/bash

g++ -std=c++11 -o3 allinone.cpp -o $1
./$1
rm ./$1