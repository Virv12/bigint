#include <cstdio>
#include "bigint.hpp"

void test_impl(size_t times, bool (*foo)(), char const str[]) {
	while (times--)
		if (!foo()) {
			puts("test failed:");
			puts(str);
			exit(1);
		}
}

#define TEST(times, ...) test_impl(times, __VA_ARGS__, #__VA_ARGS__)

int main() {
	TEST(1, [] { return std::is_trivial_v<bigint<64>>; });

	TEST(1e6, [] {
		auto const a = bigint<8>::rand();
		return a == a;
	});

	TEST(1e6, [] {
		auto const a = bigint<8>::rand();
		return a == - -a;
	});

	TEST(1e6, [] {
		auto const a = bigint<8>::rand();
		return a - a == 0;
	});

	TEST(1e6, [] {
		auto const a = bigint<8>::rand();
		auto const b = bigint<8>::rand();
		return a - b + b == a;
	});

	TEST(1e6, [] {
		auto const a = bigint<8>::rand();
		auto const b = bigint<8>::rand();
		auto const [q, r] = div(a, b);
		return r < b && q*b + r == a;
	});

	puts("All tests passed");
}
