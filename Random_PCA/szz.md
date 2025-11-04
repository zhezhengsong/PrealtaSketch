### To run RandomPCA.cpp
/opt/homebrew/bin/g++-15 -std=gnu++17 -O2 -Wall -Wextra \
  -I/opt/homebrew/include/eigen3 \
  -I"$(brew --prefix igraph)/include" \
  RandomPCA.cpp \
  -L"$(brew --prefix igraph)/lib" -ligraph \
  -o RandomPCA

./RandomPCA