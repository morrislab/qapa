#!/bin/bash

for i in python/*_test.py
do
    python $i
    echo -e "\n"
done

