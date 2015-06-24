#!/bin/bash

git submodule update --init
make -C corr
sudo make -C corr install

make
./setup.py build
sudo ./setup.py install

cd squmfit
sudo ./setup.py install
