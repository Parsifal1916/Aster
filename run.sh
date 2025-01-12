rm -rf build 
mkdir build 

cmake -B build -DCMAKE_INSTALL_PREFIX=/usr/local 
cmake --build build
sudo cmake --install build


g++ -o tests/demo tests/demo.cpp -lsfml-graphics -lsfml-window -lsfml-system -lAster g++ demo.cpp -o demo  -mtune=native -march=native -funroll-loops -ffast-math -fcx-limited-range -flto -lsfml-graphics -lsfml-window -lsfml-system -pthread -Wno-narrowing -Wno-aggressive-loop-optimizations -O3  && tests/demo
