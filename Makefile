all: run_test run_benchmark

run_test run_benchmark: run_% : %
	./$<

test: test.cpp bigint*.hpp
	clang++ -std=c++20 -Og -Wall -Wextra -Wpedantic -Wshadow -D_GLIBCXX_DEBUG -fsanitize=address -o $@ $<

benchmark: benchmark.cpp bigint*.hpp
	clang++ -std=c++20 -Ofast -Wall -Wextra -Wpedantic -Wshadow -o $@ $<

clean:
	rm -rf test benchmark

vim_1: run_test
vim_2: run_benchmark

.PHONY: all run_test run_benchmark vim_1 vim_2
