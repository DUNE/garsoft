#/bin/bash

if [ ! -d ${PWD}/build ]; then
  mkdir -p ${PWD}/build
fi

cd ${PWD}/build
rm -rf *
cmake ..
echo "done with cmake"
make -j4 install
echo "done with make install"
cd ..
