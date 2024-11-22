#pragma once
#include <vector>

// This class provides an iterator that decomposes a mode occupation multi-index into its constituents.

class occupation_numbers {
public:
	template<typename ModeIterator>
	class iterator {
	public:
		int operator*() const {
			return idx_%*mode_dim_it_;
		}
		void operator++() {
			idx_ /= *mode_dim_it_;
			mode_dim_it_++;
		}
		bool operator!=(const iterator &other) const {
			return mode_dim_it_ != other.mode_dim_it_;
		}
	private:
		iterator(int idx, ModeIterator mode_dim_it) : idx_{idx}, mode_dim_it_{mode_dim_it} {
		}
			
		int idx_{};
		ModeIterator mode_dim_it_;

		friend class occupation_numbers;
	};		
	
	auto begin() const {
		return iterator{idx_, mode_dims_.begin()};
	}

	auto end() const {
		return iterator{idx_, mode_dims_.end()};
	}
		
	occupation_numbers(int idx, const std::vector<int> &mode_dims) : idx_{idx}, mode_dims_(mode_dims) {
	}

private:
	int idx_{};
	const std::vector<int> &mode_dims_;
};
