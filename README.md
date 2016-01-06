# RegAlloc
Chaitin-Briggs register-allocation algorithm (LLVM back-end)

The Chaitin-Briggs register allocation algorithm is an interesting way of using graph coloring techniques to determine if register spilling will be necessary for a given program, and if so, which register(s) to spill. We implement this algorithm in the LLVM backend and compare the performance with LLVM-provided register allocation algorithms.
