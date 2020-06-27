#include <cstdio>
#include <vector>
#include <tuple>
#include <thread>
#include "bigint.hpp"

using namespace std::literals;

template <int> int unused;

std::vector<std::tuple<size_t, void(*)(), char const*>> foos;

#define BENCHMARK(times, foo) template<> static auto unused<__LINE__> = foos.emplace_back(times, foo, #foo)

int main() {
	for (auto &[times, foo, str] : foos) {
		std::this_thread::sleep_for(50ms);

		timespec t0, t1;
		clock_gettime(0, &t0);
		for (size_t i = 0; i < times; i++)
			foo();
		clock_gettime(0, &t1);
		printf("%gns\t%s\n", double((t1.tv_sec - t0.tv_sec) * long(1e9) + (t1.tv_nsec - t0.tv_nsec)) / times, str);
	}
}

bigint<64> a = bigint<64>::rand();
bigint<64> b = bigint<64>::rand();
bigint<64> c, d;

BENCHMARK(1e6, [] { std::tie(c, d) = div(a, b); });
BENCHMARK(1e6, [] { c = a % b; });
BENCHMARK(1e6, [] { c = a / b; });
