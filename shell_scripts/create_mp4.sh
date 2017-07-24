#! /bin/bash

ffmpeg -framerate 10 -i /tmp/output_%03.png -codec copy /tmp/test.mp4
