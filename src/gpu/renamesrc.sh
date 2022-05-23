#!/bin/bash

for src in ./*
do
src=${src##*/}
ext=${src##*.}
if [ "$ext" == "h" ] || [ "$ext" == "cpp" ] || [ "$ext" == "cuh" ] || [ "$ext" == "cu" ]; then
	mv "$src" "$(echo "$src" | sed s/pf//)"
fi
done