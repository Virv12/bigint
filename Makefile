run: main
	./main

main: main.cpp bigint.hpp
	clang++ -std=c++20 -Ofast -Wall -Wextra -Wpedantic -Wshadow -o $@ $<

clean:
	rm -rf main

vim_1: main

.PHONY: run clean vim_1
