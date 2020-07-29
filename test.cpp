#include <cstdio>
#include <iostream>
#include "bigint.hpp"

static std::mt19937_64 random_engine([] {
	auto const seed = time(0);
	std::cerr << "Test seed: " << seed << '\n';
	return seed;
} ());

template <class... T>
void test_impl(size_t times, bool (*foo)(T...), char const str[]) {
	auto test = [foo, str]<size_t... Idx>(std::index_sequence<Idx...>) {
		auto x = std::tuple{T::rand(random_engine)...};
		bool res = foo(std::get<Idx>(x)...);
		if (!res) {
			std::cerr << "Test failed: " << str << '\n';
			((std::cerr << std::get<Idx>(x) << '\n'), ...);
			exit(1);
		}
	};

	while (times--)
		test(std::make_index_sequence<sizeof...(T)>());

	std::cerr << "Test passed: " << str << '\n';
}

#define TEST(times, ...) test_impl(times, +__VA_ARGS__, #__VA_ARGS__)

int main() {
	TEST(1,   [] { return std::is_trivial_v<bigint<8>>; });
	TEST(1e6, [](bigint<8> a) { return a == a; });
	TEST(1e6, [](bigint<8> a) { return a == - -a; });

	TEST(1e6, [](bigint<8> a)                           { return a + 0 == a; });
	TEST(1e6, [](bigint<8> a, bigint<8> b)              { return a + b == b + a; });
	TEST(1e6, [](bigint<8> a, bigint<8> b, bigint<8> c) { return (a + b) + c == a + (b + c); });

	TEST(1e6, [](bigint<8> a) { return a - a == 0; });
	TEST(1e6, [](bigint<8> a, bigint<8> b) { return a - b + b == a; });

	TEST(1e6, [](bigint<8> a, bigint<8> b) {
		auto const [q, r] = div(a, b);
		return r < b && q*b + r == a;
	});
}
