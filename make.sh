#!/bin/sh

g++ --std=c++11 -O3 rtree.cpp -o rtree -I local/include/spatialindex/ -L local/lib/ -lpthread -lspatialindex
