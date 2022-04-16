#! /bin/sh

wget https://boostorg.jfrog.io/artifactory/main/release/1.75.0/source/boost_1_75_0.tar.gz -P "$(pwd)"
tar -xvf "$(pwd)"/boost_1_75_0.tar.gz

tar -xvf sdsl-lite.tar.gz
cd sdsl-lite
./install.sh "$(pwd)"/libsdsl
mv libsdsl/ ..

