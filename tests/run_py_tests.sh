#!/bin/bash

for i in python/test_*.py
do
    python $i
    echo -e "\n"
done

