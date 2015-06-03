#!/bin/bash
run:
	./my_convolution 3 3 1 -1 -2 -1 0 0 0 1 2 1

install: 
	terra convolve.lua

sobel:
	terra convolve.t && ./my_convolution 3 3 1 -1 -2 -1 0 0 0 1 2 1
		#use: output + 0.5f for each value

blur:
	terra convolve.t && ./my_convolution 5 5 273 1 4 7 4 1  4 16 26 16 4  7 26 41 26 7  4 16 26 16 4  1 4 7 4 1 
