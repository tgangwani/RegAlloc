#!/bin/bash

# makeTests.sh - Builds and optionally runs all test cases for unit project
#   ./makeTests.sh = Build all test cases
#   ./makeTests.sh run = Build all test cases, then run them (use keyboard to advance)

# To change the allocator used, change the "unitP" on line 50 to any of the other
# allocators, such as "basic", "greedy", "pbqp", "fast -O0", or "spillAll".
# Additionally, don't forget to also change the allocator for Olden/tsp on lines 58-61.

all_c_files="
./SingleSource/Benchmarks/BenchmarkGame/fannkuch
./SingleSource/Benchmarks/Misc/ReedSolomon
./SingleSource/Benchmarks/Shootout/lists
./MultiSource/Benchmarks/llubenchmark/llubenchmark
./MultiSource/Benchmarks/NPB-serial/is/is
./MultiSource/Benchmarks/Trimaran/netbench-url/packet
./MultiSource/Benchmarks/Trimaran/netbench-url/search
./MultiSource/Benchmarks/Trimaran/netbench-url/url
./MultiSource/Benchmarks/Trimaran/netbench-url/utils
./MultiSource/Benchmarks/Trimaran/enc-rc4/rc4
./MultiSource/Benchmarks/Olden/health/args
./MultiSource/Benchmarks/Olden/health/health
./MultiSource/Benchmarks/Olden/health/list
./MultiSource/Benchmarks/Olden/health/poisson
./MultiSource/Benchmarks/VersaBench/bmm/bmm
"

all_standalone_exes="
./SingleSource/Benchmarks/BenchmarkGame/fannkuch
./SingleSource/Benchmarks/Misc/ReedSolomon
./SingleSource/Benchmarks/Shootout/lists
./MultiSource/Benchmarks/llubenchmark/llubenchmark
./MultiSource/Benchmarks/NPB-serial/is/is
./MultiSource/Benchmarks/Trimaran/enc-rc4/rc4
./MultiSource/Benchmarks/VersaBench/bmm/bmm
"

no_input_exe_commands="
./SingleSource/Benchmarks/BenchmarkGame/fannkuch
./SingleSource/Benchmarks/Misc/ReedSolomon
./SingleSource/Benchmarks/Shootout/lists
./MultiSource/Benchmarks/NPB-serial/is/is
"

# Compile and run llc on all test cases:
for testfile in $all_c_files
do
	clang -O2 -emit-llvm -c $testfile.c -o $testfile.bc
	llc -march=x86 -stats -debug -regalloc=unitP $testfile.bc -o $testfile.s &> $testfile.txt
done

# The Olden/tsp test case requires special parameters to clang and llc to compile:
clang -O2 -emit-llvm -DTORONTO -c ./MultiSource/Benchmarks/Olden/tsp/args.c -o ./MultiSource/Benchmarks/Olden/tsp/args.bc
clang -O2 -emit-llvm -DTORONTO -c ./MultiSource/Benchmarks/Olden/tsp/build.c -o ./MultiSource/Benchmarks/Olden/tsp/build.bc
clang -O2 -emit-llvm -DTORONTO -c ./MultiSource/Benchmarks/Olden/tsp/main.c -o ./MultiSource/Benchmarks/Olden/tsp/main.bc
clang -O2 -emit-llvm -DTORONTO -c ./MultiSource/Benchmarks/Olden/tsp/tsp.c -o ./MultiSource/Benchmarks/Olden/tsp/tsp.bc
llc -march=x86 -stats -debug -regalloc=unitP ./MultiSource/Benchmarks/Olden/tsp/args.bc -o ./MultiSource/Benchmarks/Olden/tsp/args.s &> ./MultiSource/Benchmarks/Olden/tsp/args.txt
llc -march=x86 -stats -debug -regalloc=unitP ./MultiSource/Benchmarks/Olden/tsp/build.bc -o ./MultiSource/Benchmarks/Olden/tsp/build.s &> ./MultiSource/Benchmarks/Olden/tsp/build.txt
llc -march=x86 -stats -debug -regalloc=unitP ./MultiSource/Benchmarks/Olden/tsp/main.bc -o ./MultiSource/Benchmarks/Olden/tsp/main.s &> ./MultiSource/Benchmarks/Olden/tsp/main.txt
llc -march=x86 -stats -debug -regalloc=unitP ./MultiSource/Benchmarks/Olden/tsp/tsp.bc -o ./MultiSource/Benchmarks/Olden/tsp/tsp.s &> ./MultiSource/Benchmarks/Olden/tsp/tsp.txt

# Create exes for all test cases:
for testfile in $all_standalone_exes
do
	gcc -m32 -o $testfile.exe $testfile.s
done

# These test cases are made up of multiple files and need to be sent to gcc all at once:
gcc -m32 -o ./MultiSource/Benchmarks/Trimaran/netbench-url/url.exe ./MultiSource/Benchmarks/Trimaran/netbench-url/url.s ./MultiSource/Benchmarks/Trimaran/netbench-url/utils.s ./MultiSource/Benchmarks/Trimaran/netbench-url/packet.s ./MultiSource/Benchmarks/Trimaran/netbench-url/search.s
gcc -m32 -o ./MultiSource/Benchmarks/Olden/health/health.exe ./MultiSource/Benchmarks/Olden/health/health.s ./MultiSource/Benchmarks/Olden/health/args.s ./MultiSource/Benchmarks/Olden/health/list.s ./MultiSource/Benchmarks/Olden/health/poisson.s
gcc -m32 -o ./MultiSource/Benchmarks/Olden/tsp/main.exe ./MultiSource/Benchmarks/Olden/tsp/main.s ./MultiSource/Benchmarks/Olden/tsp/args.s ./MultiSource/Benchmarks/Olden/tsp/build.s ./MultiSource/Benchmarks/Olden/tsp/tsp.s -lm

if [ $# -eq 1 ]
  then
	
	# Run all exes that don't require parameters:
	for exefile in $no_input_exe_commands
	do
		read -p "Next file: $exefile" -n1 -s
		time $exefile.exe
	done

	# Exes that require input parameters to run correctly:
	read -p "Next file: llubenchmark" -n1 -s
	time ./MultiSource/Benchmarks/llubenchmark/llubenchmark.exe -i 3000
	read -p "Next file: url" -n1 -s
	time ./MultiSource/Benchmarks/Trimaran/netbench-url/url.exe ./MultiSource/Benchmarks/Trimaran/netbench-url/medium_inputs 5000
	read -p "Next file: rc4" -n1 -s
	time ./MultiSource/Benchmarks/Trimaran/enc-rc4/rc4.exe 2000000
	read -p "Next file: health" -n1 -s
	time ./MultiSource/Benchmarks/Olden/health/health.exe 8 500 2
	read -p "Next file: tsp" -n1 -s
	time ./MultiSource/Benchmarks/Olden/tsp/main.exe 10240000
	read -p "Next file: bmm" -n1 -s
	time ./MultiSource/Benchmarks/VersaBench/bmm/bmm.exe 1000 1000

fi

