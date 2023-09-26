#!/bin/bash

sudo apptainer build --force pst-classifier.sif pst-classifier.def

apptainer exec --bind $PWD/artifacts/:/artifacts pst-classifier.sif sh -c 'cp /PstClassifierSeqan/build/src/calculate-distances /PstClassifierSeqan/build/src/pst-score-sequences /PstClassifierSeqan/build/src/pst-classifier /PstClassifierSeqan/build/src/pst-batch-training /PstClassifierSeqan/build/src/bic /artifacts'
