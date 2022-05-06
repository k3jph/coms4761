#!/bin/bash

cd software/pytorch_DGCNN
cd lib
make -j4
cd "$(dirname "$0")"
pip3 install --user numpy
pip3 install --user scipy
pip3 install --user networkx
pip3 install --user tqdm
pip3 install --user sklearn
pip3 install --user gensim
