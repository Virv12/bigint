#ifndef BIGINT_SUM
#define BIGINT_SUM 1

#include "bigint_base.hpp"

template <size_t N>
constexpr bigint<N> &operator+= (bigint<N> &a, bigint<N> const &b) noexcept {
	uint64_t carry = 0;
	for (size_t i = 0; i < N; i++) {
		auto s = (__uint128_t)a[i] + b[i] + carry;
		a[i] = s;
		carry = s >> 64;
	}
	return a;
}

template <size_t N>
[[nodiscard]] constexpr bigint<N> &&operator+ (bigint<N> &&a, bigint<N> const &b) noexcept { return std::move(a += b); }
template <size_t N>
[[nodiscard]] constexpr bigint<N> &&operator+ (bigint<N> const &a, bigint<N> &&b) noexcept { return std::move(b += a); }
template <size_t N>
[[nodiscard]] constexpr bigint<N> &&operator+ (bigint<N> &&a, bigint<N> &&b) noexcept { return std::move(a += b); }

template <size_t N>
[[nodiscard]] constexpr bigint<N> operator+ (bigint<N> const &a, bigint<N> const &b) noexcept {
	bigint<N> r;

	uint64_t carry = 0;
	for (size_t i = 0; i < N; i++) {
		auto s = (__uint128_t)a[i] + b[i] + carry;
		r[i] = s;
		carry = s >> 64;
	}

	return r;
}

template <size_t N>
[[nodiscard]] constexpr bigint<N> operator+ (bigint<N> const &a, uint64_t carry) noexcept {
	bigint<N> r;

	for (size_t i = 0; carry && i < N; i++) {
		auto s = (__uint128_t)a[i] + carry;
		r[i] = s;
		carry = s >> 64;
	}

	return r;
}

template <size_t N>
[[nodiscard]] constexpr bigint<N> &&operator+ (bigint<N> &&a, uint64_t carry) noexcept {
	for (size_t i = 0; carry && i < N; i++) {
		auto s = (__uint128_t)a[i] + carry;
		a[i] = s;
		carry = s >> 64;
	}

	return std::move(a);
}

template <size_t N>
[[nodiscard]] constexpr bigint<N> operator+ (uint64_t a, bigint<N> const &b) noexcept { return b + a; }
template <size_t N>
[[nodiscard]] constexpr bigint<N> &&operator+ (uint64_t a, bigint<N> &&b) noexcept { return std::move(b) + a; }

template <size_t N>
[[nodiscard]] constexpr bigint<N> operator- (bigint<N> const &a) noexcept { return ~a + 1; }
template <size_t N>
[[nodiscard]] constexpr bigint<N> &&operator- (bigint<N> &&a) noexcept { return ~std::move(a) + 1; }

template <size_t N>
constexpr bigint<N> &operator-= (bigint<N> &a, bigint<N> const &b) noexcept {
	uint64_t carry = 0;
	for (size_t i = 0; i < N; i++) {
		auto t = (__int128)a[i] - b[i] - carry;
		a[i] = t;
		carry = t < 0;
	}
	return a;
}

template <size_t N>
[[nodiscard]] constexpr bigint<N> &&operator- (bigint<N> &&a, bigint<N> const &b) noexcept { return std::move(a -= b); }
template <size_t N>
[[nodiscard]] constexpr bigint<N> operator- (bigint<N> const &a, bigint<N> const &b) noexcept {
	bigint<N> r;
	uint64_t carry = 0;
	for (size_t i = 0; i < N; i++) {
		auto t = (__int128)a[i] - b[i] - carry;
		r[i] = t;
		carry = t < 0;
	}
	return r;
}

#endif
