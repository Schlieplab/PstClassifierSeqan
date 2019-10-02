#!/bin/bash

docker run -ti \
	-v $PWD/CP007136.1.fa:/data/sequences/CP007136.1.fa \
	pst-classifier-seqan:latest pst-classifier /data/sequences/CP007136.1.fa 
