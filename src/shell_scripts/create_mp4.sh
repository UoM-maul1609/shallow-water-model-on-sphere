#! /bin/bash

ffmpeg -framerate 10 -i ../../tests/output_%03d.png -codec copy ../../tests/test.mp4
