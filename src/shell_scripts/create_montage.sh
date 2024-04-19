#! /bin/bash

#for i in ../../tests/output_*.png; do mv $i ${i:0:15}${i:21:10}; done


for ((i=1;i<136;i+=1))
do
var=$(printf "%03d" $i)
montage ../../tests/output_00_$var.png ../../tests/output_01_$var.png \
	../../tests/output_02_$var.png /tmp/pics/output_03_$var.png ../../tests/output_04_$var.png \
	../../tests/output_05_$var.png -geometry +0+0 ../../tests/output_$var.png

done

#ffmpeg -framerate 10 -i output_%03.png -codec copy test.mp4
