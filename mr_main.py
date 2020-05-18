#!/usr/bin/env python
#mr_main.py
# stands for my really main file 
# where we are going to be importing and trying to make a multithreading/ multiprocessing / gpu accellerated big mess of fun this is going to be the test file for all the functions in which I Have to call and combine. 
import rscu as rs
import fileparser as fpa
import DNA_calculators as dna
import sys
import multiproccessing as mpi
import multithreading as mti





def main():
	file0=sys.argv[1]
	file1=sys.argv[2]
	file2=sys.argv[3]
	reader=fpa.file_reader(file0)
	readed=fpa.file_reader(file1)
	reads=reader.file_readers()
	reades=readed.file_readers()

	




