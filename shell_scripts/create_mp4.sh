#! /bin/bash

ffmpeg -framerate 10 -i /tmp/pics/output_%03.png -codec copy /tmp/test.mp4
