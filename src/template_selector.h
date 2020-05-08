#pragma once

// This is a questionable template trick to allow switching a couple of template parameters using
// runtime values. eldritch knowledge taken from
// https://stackoverflow.com/questions/41723704/how-to-filter-a-stdinteger-sequence
#include <tuple>

namespace template_selector_impl {

template<auto>
struct value {};

template<auto...>
struct value_sequence {};

template<auto... As, auto... Bs>
constexpr value_sequence<As..., Bs...> operator+(value_sequence<As...>, value_sequence<Bs...>) {
	return {};
}

template<auto Val, size_t R>
constexpr auto filter_single(value<Val>, value<R>) {
	if constexpr((R & (1 << Val)) != 0) {
		return value_sequence<Val>{};
	} else {
		return value_sequence<>{};
	}
}

template<auto... Vals, size_t R>
constexpr auto filter(std::index_sequence<Vals...>, value<R>) {
	return (filter_single(value<Vals>{}, value<R>{}) + ...);
}

template<typename F, typename Vals, auto... Id>
auto apply_tuple(F &func, Vals &tup, value_sequence<Id...>) {
	return func(std::get<Id>(tup)...);
}

template<typename F, typename T, size_t N = std::tuple_size_v<T>, size_t R = (1 << N) - 1>
struct checker {
	static auto apply(const F &func, T &tup, const std::array<bool, N> &args) {
		bool success = true;
		for(size_t i = 0; i < N; i++) {
			if(args[i] ^ (!!(R & (1 << i)))) {
				success = false;
				break;
			}
		}
		if(success) {
			return apply_tuple(func, tup, filter(std::make_index_sequence<N>{}, value<R>{}));
		} else {
			return checker<F, T, N, R - 1>::apply(func, tup, args);
		}
	}
};

template<typename F, typename T, size_t N>
struct checker<F, T, N, 0> {
	static auto apply(const F &func, T &tup, const std::array<bool, N> &) {
		using ret_type = decltype(apply_tuple(func, tup, value_sequence<0>{}));
		if constexpr(std::is_same<ret_type, void>::value) {
		} else {
			return ret_type{};
		}
	};
};

}

template<typename F, typename T, size_t N>
auto template_select(const F &func, T &tuple, const std::array<bool, N> &args) {
	return template_selector_impl::checker<F, T>::apply(func, tuple, args);
}
