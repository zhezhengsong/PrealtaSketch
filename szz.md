### To find best params for CountSketch_v2.cpp
python3 sweep_CountSketch.py
This works before I added time to CountSketch_v2.

### To run CountSketch_v2.cpp on my Mac
/opt/homebrew/bin/g++-15 -std=gnu++17 -O2 -Wall -Wextra -I/opt/homebrew/include -L/opt/homebrew/lib CountSketch_v2.cpp -o CountSketch_v2 -ligraph
./CountSketch_v2 1000 20 1.2

### To run CountSketch.cpp on my Mac
/opt/homebrew/bin/g++-15 -std=gnu++17 -O2 -Wall -Wextra -I/opt/homebrew/include -L/opt/homebrew/lib CountSketch.cpp -o CountSketch -ligraph
./CountSketch

### To run Eigen?
clang++ -std=c++17 -O3 try_RP.cpp -I /opt/homebrew/include/eigen3 -o try_RP
