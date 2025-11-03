# Everything here was running outside this folder before!

### To run RP_v2_time.cpp
/opt/homebrew/bin/g++-15 -std=gnu++17 -O2 -Wall -Wextra \
  -I/opt/homebrew/include/eigen3 \
  -I"$(brew --prefix igraph)/include" \
  RP_v2_time.cpp \
  -L"$(brew --prefix igraph)/lib" -ligraph \
  -o RP_v2_time

### To run RP_v2.cpp
/opt/homebrew/bin/g++-15 -std=gnu++17 -O2 -Wall -Wextra \
  -I/opt/homebrew/include/eigen3 \
  -I"$(brew --prefix igraph)/include" \
  RP_v2.cpp \
  -L"$(brew --prefix igraph)/lib" -ligraph \
  -o RP_v2
  
./RP_v2 1024 15 1.0 3

### To run RP.cpp
/opt/homebrew/bin/g++-15 -std=gnu++17 -O2 -Wall -Wextra \
  -I/opt/homebrew/include/eigen3 \
  -I"$(brew --prefix igraph)/include" \
  RP.cpp \
  -L"$(brew --prefix igraph)/lib" -ligraph \
  -o RP

./RP

### To run RP.cpp
/opt/homebrew/bin/g++-15 -std=gnu++17 -O2 -Wall -Wextra \
  -I/opt/homebrew/include/eigen3 \
  -I"$(brew --prefix igraph)/include" \
  try_RP_double.cpp \
  -L"$(brew --prefix igraph)/lib" -ligraph \
  -o try_RP_double
