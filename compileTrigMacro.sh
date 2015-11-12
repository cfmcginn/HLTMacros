#!/bin/bash

mkdir -p pdfDir
g++ $1 $(root-config --cflags --libs) -std=c++11  -Wall -O2 -o "${1/%.C/}.exe" 