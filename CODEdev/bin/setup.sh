mkdir build
cd ./build
cmake ../../src
make
cp cfd.x ../
cd ../

cp ../src/io/input/input*.dat ./
