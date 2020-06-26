#include <cstdint>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <array>
#include <ostream>
#include <bit>

template <size_t N>
struct bigint {
	[[nodiscard]] constexpr bigint() noexcept = default;
	[[nodiscard]] constexpr bigint(uint64_t v) noexcept : digits{v} {}

	template <size_t N2>
	[[nodiscard]] constexpr bigint(bigint<N2> const& oth) noexcept {
		std::copy(&oth[0], &oth[std::min(N, N2)], digits.begin());
		std::fill(digits.begin() + std::min(N2, N), digits.end(), 0);
	}

	[[nodiscard]] constexpr uint64_t& operator[] (size_t i) noexcept { return digits[i]; }
	[[nodiscard]] constexpr uint64_t const& operator[] (size_t i) const noexcept { return digits[i]; }

	[[nodiscard]] constexpr bigint operator~ () const& noexcept {
		bigint r;
		for (size_t i = 0; i < N; i++)
			r[i] = ~(*this)[i];
		return r;
	}

	[[nodiscard]] constexpr bigint&& operator~ () && noexcept {
		for (size_t i = 0; i < N; i++)
			(*this)[i] = ~(*this)[i];
		return std::move(*this);
	}

	friend auto& operator<< (std::ostream& os, bigint const& bi) {
		if constexpr (N == 0) return os << "0x0";
		constexpr char const* s = "0123456789abcdef";

		os << "0x";

		for (auto p = (uint8_t const*)(bi.digits.data() + N) - 1; p >= (uint8_t const*)bi.digits.data(); p--) {
			os << s[*p >> 4] << s[*p & 0xf];
			// if ((p - (uint8_t const*)bi.digits.data()) % 8 == 0) os << ' ';
		}

		return os;
	}

	constexpr bigint& operator<<= (size_t n) noexcept {
		size_t const upper = n / 64;
		size_t const lower = n % 64;

		if (!lower) {
			std::copy(digits.rbegin() + upper, digits.rend(), digits.rbegin());
		} else {
			for (size_t i = N - 1; i-upper-1 < N; i--)
				digits[i] = (digits[i-upper] << lower) | (digits[i-upper-1] >> (64 - lower));
			digits[upper] = digits[0] << lower;
		}
		std::fill(digits.begin(), digits.begin() + upper, 0);

		return *this;
	}

	[[nodiscard]] constexpr bigint&& operator<< (size_t n) && noexcept { return std::move(*this <<= n); }
	[[nodiscard]] constexpr bigint operator<< (size_t n) const& noexcept {
		size_t const upper = n / 64;
		size_t const lower = n % 64;
		bigint r;

		if (!lower) {
			std::copy(digits.rbegin() + upper, digits.rend(), r.digits.rbegin());
		} else {
			for (size_t i = N - 1; i-upper-1 < N; i--)
				r.digits[i] = (digits[i-upper] << lower) | (digits[i-upper-1] >> (64 - lower));
			r.digits[upper] = digits[0] << lower;
		}
		std::fill(r.digits.begin(), r.digits.begin() + upper, 0);

		return r;
	}

	constexpr bigint& operator>>= (size_t n) noexcept {
		size_t const upper = n / 64;
		size_t const lower = n % 64;

		if (!lower) {
			std::copy(digits.begin() + upper, digits.end(), digits.begin());
		} else {
			for (size_t i = 0; i + upper + 1 < N; ++i)
				digits[i] = (digits[i + upper] >> lower) | (digits[i + upper + 1] << (64 - lower));
			digits[N - upper - 1] = digits[N-1] >> lower;
		}
		std::fill(digits.begin() + N - upper, digits.end(), 0);

		return *this;
	}

	[[nodiscard]] constexpr bigint&& operator>> (size_t n) && noexcept { return std::move(*this >>= n); }
	[[nodiscard]] constexpr bigint operator>> (size_t n) const& noexcept {
		size_t const upper = n / 64;
		size_t const lower = n % 64;
		bigint r;

		if (!lower) {
			std::copy(digits.begin() + upper, digits.end(), r.digits.begin());
		} else {
			for (size_t i = 0; i + upper + 1 < N; ++i)
				r.digits[i] = (digits[i + upper] >> lower) | (digits[i + upper + 1] << (64 - lower));
			r.digits[N - upper - 1] = digits[N-1] >> lower;
		}
		std::fill(r.digits.begin() + N - upper, r.digits.end(), 0);

		return r;
	}

	[[nodiscard]] constexpr bigint operator- () const& noexcept { return ~*this + 1; }
	[[nodiscard]] constexpr bigint&& operator- () && noexcept { return ~std::move(*this) + 1; }

	constexpr bigint& operator+= (bigint const& oth) noexcept {
		uint64_t carry = 0;
		for (size_t i = 0; i < N; i++) {
			auto s = (__uint128_t)(*this)[i] + oth[i] + carry;
			(*this)[i] = s;
			carry = s >> 64;
		}
		return *this;
	}

	[[nodiscard]] friend constexpr bigint&& operator+ (bigint&& a, bigint const& b) noexcept { return std::move(a += b); }
	[[nodiscard]] friend constexpr bigint&& operator+ (bigint const& a, bigint&& b) noexcept { return std::move(b += a); }
	[[nodiscard]] friend constexpr bigint&& operator+ (bigint&& a, bigint&& b) noexcept { return std::move(a += b); }
	[[nodiscard]] friend constexpr bigint operator+ (bigint const& a, bigint const& b) noexcept {
		bigint r;

		uint64_t carry = 0;
		for (size_t i = 0; i < N; i++) {
			auto s = (__uint128_t)a[i] + b[i] + carry;
			r[i] = s;
			carry = s >> 64;
		}

		return r;
	}

	constexpr bigint& operator-= (bigint const& oth) noexcept {
		uint64_t carry = 0;
		for (size_t i = 0; i < N; i++) {
			auto t = (__int128)(*this)[i] - oth[i] - carry;
			(*this)[i] = t;
			carry = t < 0;
		}
		return *this;
	}

	[[nodiscard]] constexpr bigint&& operator- (bigint const& oth) && noexcept { return std::move(*this -= oth); }
	[[nodiscard]] friend constexpr bigint operator- (bigint const& a, bigint const& b) noexcept {
		bigint r;
		uint64_t carry = 0;
		for (size_t i = 0; i < N; i++) {
			auto t = (__int128)a[i] - b[i] - carry;
			r[i] = t;
			carry = t < 0;
		}
		return r;
	}

	[[nodiscard]] friend constexpr bigint naive_mult(bigint const& a, bigint const& b) noexcept {
		bigint r;
		__uint128_t carry = 0;

		for (size_t i = 0; i < N; i++) {
			r[i] = carry;
			carry >>= 64;

			for (size_t j = 0; j < i+1; j++) {
				auto s = (__uint128_t)a.digits[j] * b.digits[i - j] + r[i];
				r[i] = s;
				carry += s >> 64;
			}
		}
		return r;
	}

	[[nodiscard]] friend constexpr bigint<2 * N> naive_mult_ext(bigint const& a, bigint const& b) noexcept {
		bigint<2 * N> r;
		__uint128_t carry = 0;

		for (size_t i = 0; i < 2 * N; i++) {
			r[i] = carry;
			carry >>= 64;

			for (size_t j = i + 1 >= N ? i - N + 1 : 0; j < std::min(i+1, N); j++) {
				auto s = (__uint128_t)a.digits[j] * b.digits[i - j] + r[i];
				r[i] = s;
				carry += s >> 64;
			}
		}

		return r;
	}

	[[nodiscard]] friend constexpr bigint<2 * N> karatsuba_ext(bigint const& a, bigint const& b) noexcept {
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

	[[nodiscard]] friend constexpr bigint karatsuba(bigint const& a, bigint const& b) noexcept {
		constexpr size_t half = (N + 1) / 2;

		bigint<half> const& a0 = *(bigint<half> const*)&a;
		bigint<half> const& b0 = *(bigint<half> const*)&b;
		bigint<half> a1 = a >> (64 * half);
		bigint<half> b1 = b >> (64 * half);

		bigint<N> x = mult_ext(a0, b0);
		bigint<N> z = mult_ext(a0, b1) + mult_ext(a1, b0);

		return (std::move(z) << (half * 64)) + x;
	}

	[[nodiscard]] friend constexpr bigint operator* (bigint const& a, bigint const& b) noexcept {
		if constexpr (true) return naive_mult(a, b);
		else                return karatsuba(a, b);
	}

	[[nodiscard]] friend constexpr bigint<2 * N> mult_ext(bigint const& a, bigint const& b) noexcept {
		if constexpr (true) return naive_mult_ext(a, b);
		else                return karatsuba_ext(a, b);
	}

	constexpr bigint& operator*= (bigint const& oth) noexcept {
		return *this = *this * oth;
	}

	[[nodiscard]] constexpr bool operator== (bigint const& b) const noexcept {
		return !memcmp(digits.data(), b.digits.data(), sizeof(digits));
	}

	[[nodiscard]] constexpr bool operator< (bigint const& b) const noexcept {
		for (size_t i = N - 1; i < N; i--) {
			if (digits[i] < b[i]) return true;
			if (digits[i] > b[i]) return false;
		}
		return false;
	}

	[[nodiscard]] constexpr size_t leading_zeros() const noexcept {
		size_t lz = 0;
		while (lz+1 < N && !digits[N - lz - 1]) lz++;
		return lz * 64 + std::countl_zero(digits[N - lz - 1]);
	}
	
	[[nodiscard]] constexpr size_t trailing_zeros() const noexcept {
		size_t lz = 0;
		while (lz+1 < N && !digits[lz]) lz++;
		return lz * 64 + std::countr_zero(digits[lz]);
	}

	[[nodiscard]] constexpr bigint operator/ (bigint const& b) const& noexcept { return div(*this, b).first; }
	[[nodiscard]] constexpr bigint operator/ (bigint const& b) && noexcept { return div(std::move(*this), b).first; }

	[[nodiscard]] constexpr bigint operator% (bigint const& b) const& noexcept { return div(*this, b).second; }
	[[nodiscard]] constexpr bigint operator% (bigint const& b) && noexcept { return div(std::move(*this), b).second; }

	// [[nodiscard]] constexpr bigint operator% (bigint const& d) const& noexcept {
	// 	bigint r = *this;

	// 	for (size_t lz = d.leading_zeros() - r.leading_zeros(); lz < 64*N; lz--) {
	// 		bigint t = d << lz;
	// 		if (!(r < t)) r -= t;
	// 	}

	// 	return r;
	// }

	// [[nodiscard]] constexpr bigint&& operator% (bigint const& d) && noexcept {
	// 	bigint&& r = std::move(*this);

	// 	for (size_t lz = d.leading_zeros() - r.leading_zeros(); lz < 64*N; lz--) {
	// 		bigint t = d << lz;
	// 		if (!(r < t)) r -= t;
	// 	}

	// 	return std::move(r);
	// }

	// [[nodiscard]] friend constexpr std::pair<bigint, bigint> div(bigint&& r, bigint const& d) noexcept {
	// 	bigint q = 0;

	// 	for (size_t lz = d.leading_zeros() - r.leading_zeros(); lz < 64*N; lz--) {
	// 		bigint t = d << lz;
	// 		if (!(r < t)) {
	// 			r -= t;
	// 			q[lz / 64] |= 1ull << (lz % 64);
	// 		}
	// 	}

	// 	return {q, r};
	// }
	
	[[nodiscard]] constexpr bigint shift(int64_t const n) const& noexcept {
		if (n < 0)
			return *this >> -n;
		return *this << n;
	}

	[[nodiscard]] constexpr bigint shift(int64_t const n) && noexcept {
		if (n < 0)
			return std::move(*this) >> -n;
		return std::move(*this) << n;
	}

	[[nodiscard]] friend constexpr std::pair<bigint, bigint> div(bigint const& ref_r, bigint const& d) noexcept {
		if constexpr (false) {
			bigint q = 0;
			bigint r = ref_r;

			int64_t dlz = 64 * N - d.leading_zeros() - 64;
			auto dr = (__uint128_t)d.shift(-dlz)[0] + 1;

			for (size_t i = N - 2; i < N; i--)
				for (size_t j = 0; j < 2 && r[i + 1]; ++j) {
					auto s = *(__uint128_t*)&r[i] / dr;
					bigint x = s;
					x[1] = s >> 64;

					x = std::move(x).shift(64 * i - dlz);

					q += x;
					r -= d * x;
				}

			if (!(r < d)) {
				uint64_t s = r[0] / d[0];
				q += s;
				r -= s * d[0];
			}

			return {q, r};

		} else {
			bigint q = 0;
			bigint r = ref_r;

			for (size_t lz = d.leading_zeros() - r.leading_zeros(); lz < 64*N; lz--) {
				bigint t = d << lz;
				if (!(r < t)) {
					r -= t;
					q[lz / 64] |= 1ull << (lz % 64);
				}
			}

			return {q, r};
		}
	}

	[[nodiscard]] constexpr bigint modinv(bigint const& mod) const {
		bigint a = *this;
		bigint b = 1;
		bigint c = 0;
		bigint n = mod;

		while (n != 0) {
			auto [q, r] = div(std::move(a), n);
			a = n;
			n = r;
			bigint s = b - q * c;
			b = c;
			c = s;
		}

		if (a != 1)
			throw std::runtime_error("base is not invertible for the given modulus");
		if (b[N-1] >> 63 & 1)
			b += mod;
		return std::move(b);
	}

	[[nodiscard]] friend constexpr bigint pow(bigint const& ref_a, bigint const& n, bigint const& m) noexcept {
		size_t const max = 64 * N - n.leading_zeros();
		bigint r = 1;
		bigint a = ref_a;

		for (size_t i = 0; i < max; i++) {
			if (n[i / 64] >> (i % 64) & 1)
				r = mult_ext(r, a) % m;
			a = mult_ext(a, a) % m;
		}

		return r;
	}

	[[nodiscard]] friend constexpr bigint pow(bigint&& a, bigint const& n, bigint const& m) noexcept {
		size_t const max = 64 * N - n.leading_zeros();
		bigint r = 1;

		for (size_t i = 0; i < max; i++) {
			if (n[i / 64] >> (i % 64) & 1)
				r = mult_ext(r, a) % m;
			a = mult_ext(a, a) % m;
		}

		return r;
	}
	
	[[nodiscard]] static bigint rand() noexcept {
		bigint r;
		for (size_t i = 0; i < N; i++)
			r[i] = ((uint64_t)std::rand() << 31 | std::rand()) << 31 | std::rand();
		return r;
	}

	[[nodiscard]] static bigint rand(bigint const& max) noexcept {
		bigint r;
		for (size_t i = 64 * N - max.leading_zeros() - 1; i < 64 * N; i--)
			r[i / 64] |= uint64_t(std::rand() & 1) << (i % 64);
		if (!(r < max))
			r -= max;
		return r;
	}

	[[nodiscard]] static bigint rand(bigint const& min, bigint const& max) noexcept {
		return rand(max - min) + min;
	}

private:
	std::array<uint64_t, N> digits;
};
