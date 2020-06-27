#ifndef BIGINT_MULT
#define BIGINT_MULT 1

#include "bigint_base.hpp"

template <size_t N>
[[nodiscard]] constexpr bigint<N> naive_mult(bigint<N> const &a, bigint<N> const &b) noexcept {
	bigint<N> r;
	__uint128_t carry = 0;

	for (size_t i = 0; i < N; i++) {
		r[i] = carry;
		carry >>= 64;

		for (size_t j = 0; j < i+1; j++) {
			auto s = (__uint128_t)a[j] * b[i - j] + r[i];
			r[i] = s;
			carry += s >> 64;
		}
	}
	return r;
}

template <size_t N>
[[nodiscard]] constexpr bigint<2 * N> naive_mult_ext(bigint<N> const &a, bigint<N> const &b) noexcept {
	bigint<2 * N> r;
	__uint128_t carry = 0;

	for (size_t i = 0; i < 2 * N; i++) {
		r[i] = carry;
		carry >>= 64;

		for (size_t j = i + 1 >= N ? i - N + 1 : 0; j < std::min(i+1, N); j++) {
			auto s = (__uint128_t)a[j] * b[i - j] + r[i];
			r[i] = s;
			carry += s >> 64;
		}
	}

	return r;
}

template <size_t N>
[[nodiscard]] constexpr bigint<2 * N> karatsuba_ext(bigint<N> const &a, bigint<N> const &b) noexcept {
	constexpr size_t half = (N + 1) / 2;

	bigint<half+1> a0 = *(bigint<half> const*)&a;
	bigint<half+1> b0 = *(bigint<half> const*)&b;
	bigint<half+1> a1 = a >> (64 * half);
	bigint<half+1> b1 = b >> (64 * half);

	if constexpr (true) {
		bigint<2 * N> r;

		auto x = mult_ext(a0, b0);
		std::copy(&x[0], &x[N], &r[0]);

		auto y = mult_ext(b0, b1);
		std::copy(&y[0], &y[N], &r[N]);

		auto z = mult_ext(std::move(a0) + a1, std::move(b0) + b1);
		z -= x;
		z -= y;

		return std::move(r) + ((bigint<2 * N>)z << (64 * half));

	} else {
		bigint<2 * N> x = mult_ext(a0, b0);
		bigint<2 * N> y = mult_ext(a1, b1);
		bigint<2 * N> z = mult_ext(std::move(a0) + a1, std::move(b0) + b1);

		z = std::move(z) - x - y;
		return (((std::move(y) << (half * 64)) + std::move(z)) << (half * 64)) + x;
	}
}

template <size_t N>
[[nodiscard]] constexpr bigint<N> karatsuba(bigint<N> const &a, bigint<N> const &b) noexcept {
	constexpr size_t half = (N + 1) / 2;

	bigint<half> const &a0 = *(bigint<half> const*)&a;
	bigint<half> const &b0 = *(bigint<half> const*)&b;
	bigint<half> a1 = a >> (64 * half);
	bigint<half> b1 = b >> (64 * half);

	bigint<N> x = mult_ext(a0, b0);
	bigint<N> z = mult_ext(a0, b1) + mult_ext(a1, b0);

	return (std::move(z) << (half * 64)) + x;
}

template <size_t N>
[[nodiscard]] constexpr bigint<N> operator* (bigint<N> const &a, bigint<N> const &b) noexcept {
	if constexpr (true) return naive_mult(a, b);
	else                return karatsuba(a, b);
}

template <size_t N>
[[nodiscard]] constexpr bigint<2 * N> mult_ext(bigint<N> const &a, bigint<N> const &b) noexcept {
	if constexpr (true) return naive_mult_ext(a, b);
	else                return karatsuba_ext(a, b);
}

template <size_t N>
constexpr bigint<N> &operator*= (bigint<N> &a, bigint<N> const &b) noexcept {
	return a = a * b;
}

#endif
