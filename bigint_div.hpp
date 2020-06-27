#ifndef BIGINT_DIV
#define BIGINT_DIV 1

#include <utility>
#include "bigint_base.hpp"

template <size_t N>
[[nodiscard]] constexpr bigint<N> operator/ (bigint<N> const &a, bigint<N> const &b) noexcept { return div(a, b).first; }
template <size_t N>
[[nodiscard]] constexpr bigint<N> operator/ (bigint<N> &&a, bigint<N> const &b) noexcept { return div(std::move(a), b).first; }

template <size_t N>
[[nodiscard]] constexpr bigint<N> operator% (bigint<N> const &r_ref, bigint<N> const &d) noexcept {
	bigint<N> r = r_ref;

	for (size_t lz = d.leading_zeros() - r.leading_zeros(); lz < 64*N; lz--) {
		bigint<N> t = d << lz;
		if (!(r < t)) r -= t;
	}

	return r;
}

template <size_t N>
[[nodiscard]] constexpr bigint<N> &&operator% (bigint<N> &&r, bigint<N> const &d) noexcept {
	for (size_t lz = d.leading_zeros() - r.leading_zeros(); lz < 64*N; lz--) {
		bigint<N> t = d << lz;
		if (!(r < t)) r -= t;
	}

	return std::move(r);
}

template <size_t N>
[[nodiscard]] constexpr std::pair<bigint<N>, bigint<N>> div(bigint<N> const &ref_r, bigint<N> const &d) noexcept {
	bigint<N> q = 0;
	bigint<N> r = ref_r;

	for (size_t lz = d.leading_zeros() - r.leading_zeros(); lz < 64*N; lz--) {
		bigint<N> t = d << lz;
		if (!(r < t)) {
			r -= t;
			q[lz / 64] |= 1ull << (lz % 64);
		}
	}

	return {q, r};
}

template <size_t N>
[[nodiscard]] constexpr std::pair<bigint<N>, bigint<N>> div(bigint<N> &&r, bigint<N> const &d) noexcept {
	bigint<N> q = 0;

	for (size_t lz = d.leading_zeros() - r.leading_zeros(); lz < 64*N; lz--) {
		bigint<N> t = d << lz;
		if (!(r < t)) {
			r -= t;
			q[lz / 64] |= 1ull << (lz % 64);
		}
	}

	return {q, r};
}

#endif
