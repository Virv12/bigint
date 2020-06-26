#include <cstdio>
#include <thread>
#include "bigint.hpp"

using namespace std::literals;

void benchmark(size_t times, void (*foo)(), char const s[]) {
	std::this_thread::sleep_for(50ms);

	timespec t0, t1;
	clock_gettime(0, &t0);
	for (size_t i = 0; i < times; i++)
		foo();
	clock_gettime(0, &t1);
	printf("%gns\t%s\n", double((t1.tv_sec - t0.tv_sec) * long(1e9) + (t1.tv_nsec - t0.tv_nsec)) / times, s);
}

#define BENCHMARK(times, foo) benchmark(times, foo, #foo)

bigint<64> a = bigint<64>::rand();
bigint<64> b = bigint<64>::rand();
bigint<64> c;

int main() {
	BENCHMARK(1e6, [] { c = a + b; });
}
