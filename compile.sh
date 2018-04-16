# Use the first line to compile with multi-threading support (requires
# pthread and associated libraries). Use the second line to compile
# a single-thread application.
g++ -O4 -lm -fopenmp imageProc.c TimeLapseFusion.c -o TimeLapseFusion
#g++ -O4 -lm imageProc.c TimeLapseFusion.c -o TimeLapseFusion
