#ifndef BIGINT_RAND
#define BIGINT_RAND 1

#include <random>

template <size_t N>
bigint<N> bigint<N>::rand() noexcept {
	bigint<N> r;
	for (size_t i = 0; i < N; i++)
		r[i] = ((uint64_t)std::rand() << 31 | std::rand()) << 31 | std::rand();
	return r;
}

template <size_t N>
bigint<N> bigint<N>::rand(bigint<N> const &max) noexcept {
	bigint<N> r;
	for (size_t i = 64 * N - max.leading_zeros() - 1; i < 64 * N; i--)
		r[i / 64] |= uint64_t(std::rand() & 1) << (i % 64);
	if (!(r < max))
		r -= max;
	return r;
}

template <size_t N>
bigint<N> bigint<N>::rand(bigint<N> const &min, bigint<N> const &max) noexcept {
	return rand(max - min) + min;
}

template <size_t N>
bigint<N> bigint<N>::rand(std::mt19937_64 &random) {
	bigint<N> r;
	for (size_t i = 0; i < N; i++)
		r[i] = random();
	return r;
}

#endif
