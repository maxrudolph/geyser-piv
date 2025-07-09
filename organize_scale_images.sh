#!/bin/bash

for i in `find . -name img_00001.png`
do
echo $i
cp $i `echo $i | cut -d / -f 3`_scale.png
done
