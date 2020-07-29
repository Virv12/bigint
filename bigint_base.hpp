#ifndef BIGINT
#define BIGINT 1

#include <cstdint>
#include <cstring>
#include <algorithm>
#include <array>
#include <ostream>
#include <random>
#include <bit>

template <size_t N>
struct bigint {
	std::array<uint64_t, N> digits;

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

	[[nodiscard]] explicit constexpr operator bool() const noexcept {
		return std::find_if(digits.begin(), digits.end(), [] (uint64_t d) -> bool { return d; }) != digits.end();
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
				r = mult_ext(r, a) % (bigint<2*N>)m;
			a = mult_ext(a, a) % (bigint<2*N>)m;
		}

		return r;
	}

	[[nodiscard]] friend constexpr bigint pow(bigint&& a, bigint const& n, bigint const& m) noexcept {
		size_t const max = 64 * N - n.leading_zeros();
		bigint r = 1;

		for (size_t i = 0; i < max; i++) {
			if (n[i / 64] >> (i % 64) & 1)
				r = mult_ext(r, a) % (bigint<2*N>)m;
			a = mult_ext(a, a) % (bigint<2*N>)m;
		}

		return r;
	}
	
	[[deprecated("consider using bigint::rand(std::mt19937_64&)")]] [[nodiscard]] static bigint rand() noexcept;
	[[nodiscard]] static bigint rand(bigint const&) noexcept;
	[[nodiscard]] static bigint rand(bigint const&, bigint const&) noexcept;

	[[nodiscard]] static bigint rand(std::mt19937_64&);
};

#endif
