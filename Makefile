#!/bin/bash
TERRA = terra

all: sobel;

install: 
	$(TERRA) imageconv.lua

run:
	my_convolution 3 3 1 -1 -2 -1 0 0 0 1 2 1

# Examles of use
sobel:
	$(TERRA) imageconv.t && ./bin/my_imageconv 3 3 1 -1 -2 -1 0 0 0 1 2 1
		#use: output + 0.5f for each value

blur:
	$(TERRA) imageconv.t && ./bin/my_imageconv 5 5 273 1 4 7 4 1  4 16 26 16 4  7 26 41 26 7  4 16 26 16 4  1 4 7 4 1 