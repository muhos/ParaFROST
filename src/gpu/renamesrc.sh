#!/bin/bash

for src in ./*
do
file=${src##*/}
ext=${file##*.}
name=${src%.*}
if [ "$ext" == "h" ]; then
	mv "$src" "${name}.hpp"
fi
done