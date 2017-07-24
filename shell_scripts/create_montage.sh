#! /bin/bash

#for i in /tmp/pics/output_*.png; do mv $i ${i:0:15}${i:21:10}; done


for ((i=1;i<136;i+=1))
do
var=$(printf "%03d" $i)
montage /tmp/pics/output_00_$var.png /tmp/pics/output_01_$var.png \
	/tmp/pics/output_02_$var.png /tmp/pics/output_03_$var.png /tmp/pics/output_04_$var.png \
	/tmp/pics/output_05_$var.png -geometry +0+0 /tmp/pics/output_$var.png

done

#ffmpeg -framerate 10 -i output_%03.png -codec copy test.mp4
