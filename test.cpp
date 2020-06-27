#include <cstdio>
#include "bigint.hpp"

void test_impl(bool val, char const str[]) {
	if (!val) {
		puts("test failed:");
		puts(str);
		exit(1);
	}
}

#define TEST(...) test_impl(__VA_ARGS__, #__VA_ARGS__)

int main() {
	TEST(std::is_trivial_v<bigint<64>>);

	TEST([] {
		auto a = bigint<8>::rand();
		auto b = bigint<8>::rand();
		auto [q, r] = div(a, b);
		return r < b && q*b + r == a;
	} ());

	puts("All tests passed");
}
