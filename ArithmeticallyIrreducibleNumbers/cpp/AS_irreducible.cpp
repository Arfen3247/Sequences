#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <cassert>
#include <chrono>
#include <format>
#include <iostream>
#include <numeric>
#include <ranges>
#include <string>
#include <vector>

std::size_t find_upper_bound_L(const std::size_t base, bool conjecture_L) {
	if (conjecture_L){ // conjectured upper bound 2^n-1, where n=ceil(log_2(base))
		std::size_t x = 1;
		while (x < base) {
			x <<= 1;
		}
		return x - 1;
	}
	return 2 * base - 1; // safe upper bound
}

std::array<std::size_t, 3> find_basics(const std::size_t base, bool conjecture_L) {
	const auto L_upper = find_upper_bound_L(base, conjecture_L);
	const auto bound = (base - 1) * L_upper;	// upper bound on sum of irreducible string
	const auto left_bits = (bound + base - 1) / 2;
	const auto right_bits = (bound - 1) / 2;
	const auto total_bits = left_bits + 1 + right_bits;	// total bits	
	return {L_upper, right_bits, total_bits};	
}

// serial, O(b^3) memory, based on Back-Tracking
std::vector<std::size_t> count_serial_backtrack(const std::size_t base, bool conjecture_L = false) {
	auto [L_upper, right_bits, total_bits] = find_basics(base, conjecture_L);
	std::vector<std::size_t> counts(L_upper+2, 0);	// required result

	// bitmasks representing the digits not allowed to appear next
	std::vector<boost::dynamic_bitset<>> list_of_disallowed_masks;
	list_of_disallowed_masks.reserve(L_upper+2);
	list_of_disallowed_masks.emplace_back(total_bits, 0);

	// the stack for the back-tracking, with a fake initial state
	std::vector<std::size_t> stack(1,1);
	stack.reserve((base-1) * L_upper);

	// back-tracking: pop the top digit of the stack,
	// if 0, then back-track, else use this digit,
	// add possible next digits to the stack
	while (!stack.empty()) {
		std::size_t digit = stack.back();
		stack.pop_back();
		if (!digit) { 
			list_of_disallowed_masks.pop_back();
			continue;
		}
		stack.emplace_back(0);
		auto mask = list_of_disallowed_masks.back();
		mask = (mask << digit) | (mask >> digit);
		mask.set(right_bits);
		std::size_t found = 0;
		for (std::size_t i = 1; i < base; i++){
			if (!mask.test(right_bits+i)){
				stack.emplace_back(i);
				++found;
			}
		}
		counts[list_of_disallowed_masks.size()] += found;
		list_of_disallowed_masks.push_back(mask);
	}

	// remove trailing zeros
	while (!counts.back()) {
		counts.pop_back();
	}
	return counts;
}

// serial, O(b^4) memory, based on a Depth-First-Seach
// I don't see how to parallelise this without a shared memory, which sounds slow
std::vector<std::size_t> count_serial_dfs(const std::size_t base, bool conjecture_L = false) {
	auto [L_upper, right_bits, total_bits] = find_basics(base, conjecture_L);
	std::vector<std::size_t> counts(L_upper+2, 0);	// required result

	// the stack for the depth-first-search, with a fake initial state
	// each state is a bitmask marking the digits not allowed to follow, and a 'depth' counter
	typedef std::pair<boost::dynamic_bitset<>, std::size_t> stackpair;
	std::vector<stackpair> stack;
	stack.emplace_back(boost::dynamic_bitset(total_bits, 0), 0);
	stack.reserve((base-1) * L_upper);

	// DFS: pop the state on top of the stack, append its children states to the stack
	while (!stack.empty()) {
		auto [mask, depth] = std::move(stack.back());
		stack.pop_back();
		mask.set(right_bits);
		++depth;
		std::size_t found = 0;
		for (std::size_t i = 1; i < base; i++){
			if (!mask.test(right_bits+i)){
				++found;
				stack.emplace_back((mask << i) | (mask >> i), depth);
			}
		}
		counts[depth] += found;    
	}

	// remove trailing zeros
	while (!counts.back()) {
		counts.pop_back();
	}
	return counts;
}

// Parallelisation idea:
// partition the tree so no branch is relatively large,
// each branch is counted using the above Back-Tracking approach,
// separate branches can be computed in parallel

// same logic as count_serial_backtrack, but it collects the valid arrays in a vector
// should only be used for small bases, eg <= 10
std::vector<std::vector<std::size_t>> collect_serial_backtrack(const size_t base, const bool conjecture_L = false) {
	auto [L_upper, right_bits, total_bits] = find_basics(base, conjecture_L);
	std::vector<std::vector<std::size_t>> collection;	// required result
	std::vector<std::size_t> arr;	// vector representing the branch
	arr.reserve(L_upper);

	std::vector<boost::dynamic_bitset<>> list_of_disallowed_masks;
	list_of_disallowed_masks.reserve(L_upper+2);
	boost::dynamic_bitset<> mask(total_bits);
	mask.set(right_bits);
	list_of_disallowed_masks.push_back(mask);

	std::vector<std::size_t> stack;
	stack.reserve((base-1) * L_upper);
	for (std::size_t i = 1; i < base; i++) {
		stack.push_back(i);
	}

	while (!stack.empty()) {
		std::size_t digit = stack.back();
		stack.pop_back();
		if (!digit) { 
			list_of_disallowed_masks.pop_back();
			arr.pop_back();
			continue;
		}
		stack.emplace_back(0);
		arr.push_back(digit);
		collection.push_back(arr);	

		auto mask = list_of_disallowed_masks.back();
		mask = (mask << digit) | (mask >> digit);
		mask.set(right_bits);
		for (std::size_t i = 1; i < base; i++){
			if (!mask.test(right_bits + i)){
				stack.emplace_back(i);
			}
		}
		list_of_disallowed_masks.push_back(mask);
	}
	return collection;
}

// How to partition:
// choose an integer delta > 1, define scale = floor((base-1)/delta).
// If a string x is irreducible in 'base' =  scale, 
// then y = scale * x is irreducible in 'base' = base.
// Heuristic: if delta is a power of 2, the branch starting at y is expensive.
// So, all the children of y should be made separate starting points.
// The choice delta = 8 should work well for each 16 <= base <= 32
std::vector<std::vector<std::size_t>> partition_tree(const size_t base, const size_t delta = 8, const bool conjecture_L = false) {
	auto scale = (base - 1) / delta + 1;
	std::vector<std::size_t> other_digits;  // those not divisible by delta
	for (std::size_t i = 1; i < base; i++) {
		if (i % delta) {
			other_digits.push_back(i);
		}
	}	
	auto templates = collect_serial_backtrack(scale, conjecture_L);
	for (auto& t : templates) {
		for (auto& x : t) {
			x *= delta;
		}
	}
	std::vector<std::vector<std::size_t>> branches;
	branches.reserve(templates.size() * other_digits.size());
	for (const auto& t : templates) {
		for (auto& i : other_digits) {
			branches.push_back(t);
			branches.back().push_back(i);
		}
	}

	for (auto& i : other_digits) {
		std::vector<std::size_t> v;
		v.push_back(i);
		branches.push_back(v);
	}
	return branches;
}

// the original back-tracking algorithm applied to the given branch
std::vector<std::size_t> count_branch(const std::size_t base, const std::array<std::size_t, 3> basics, std::vector<std::size_t> branch) {
	const auto [L_upper, right_bits, total_bits] = basics;
	std::vector<std::size_t> counts(L_upper+2, 0);
	counts[branch.size()] = 1;

	std::vector<boost::dynamic_bitset<>> list_of_disallowed_masks;
	list_of_disallowed_masks.reserve(L_upper+2);
	boost::dynamic_bitset mask(total_bits, 0);
	list_of_disallowed_masks.push_back(mask);
	mask.set(right_bits);
	list_of_disallowed_masks.push_back(mask);

	auto prev_digit = branch.back();
	branch.pop_back();
	for (auto & digit : branch) {
		auto mask = list_of_disallowed_masks.back();
		mask = (mask << digit) | (mask >> digit);
		mask.set(right_bits);
		list_of_disallowed_masks.push_back(mask);
	}

	std::vector<std::size_t> stack(1, prev_digit);
	stack.reserve((base-1) * (L_upper - branch.size()));	//??
	while (!stack.empty()) {
		std::size_t digit = stack.back();
		stack.pop_back();
		if (!digit) { 
			list_of_disallowed_masks.pop_back();
			continue;
		}
		stack.emplace_back(0);
		auto mask = list_of_disallowed_masks.back();
		mask = (mask << digit) | (mask >> digit);
		mask.set(right_bits);
		std::size_t found = 0;
		for (std::size_t i = 1; i < base; i++){
			if (!mask.test(right_bits+i)){
				stack.emplace_back(i);
				++found;
			}
		}
		counts[list_of_disallowed_masks.size()] += found;
		list_of_disallowed_masks.push_back(mask);
	}
	return counts;
}

// parallelise this!
std::vector<std::size_t> count_parallel_backtrack(const size_t base, const std::size_t delta = 4, const bool conjecture_L = false) {
	if (base < 10) {
		return count_serial_backtrack(base, conjecture_L);	// fast enough, so we focus on parallelising 'big' bases
	}

	auto branches = partition_tree(base, delta, conjecture_L);
	auto basics = find_basics(base, conjecture_L);
	std::vector<std::size_t> counts(basics[0]+2, 0);

	auto new_counts = count_serial_backtrack((base-1)/delta+1, conjecture_L); // strange fix
	for (std::size_t i = 0; i < new_counts.size(); i++) {
		counts[i] += new_counts[i];
	}

	for (auto branch : branches) {
		auto new_counts = count_branch(base, basics, branch); 	
		for (std::size_t i = 0; i < counts.size(); i++) {
			counts[i] += new_counts[i];
		}
	}	
	while (!counts.back()) {
		counts.pop_back();
	}
	return counts;
}

int main(int argc, char *argv[]) {

	std::size_t begin, end;

	if (argc == 2 && (std::string(argv[1]) == "help" || std::string(argv[1]) == "h")) {
		std::cout << "This binary is designed to be used in 3 different modes. All argumets to this binary call are bases "
			"to be computed by the counting function.\n"
			"   1. When ran with one base argument, that base is computed.\n"
			"   2. When ran with two base arguments an inclusive range of those two bases is computed.\n"
			"   3. When ran with no arguments a default range of 2 to 10 is computed.\n"
			"\n"
			"Examples on how to run this binary\n"
			"\n"
			"    ./AS_reducible 5\n"
			"    ./AS_reducible 3 10\n"
			"    ./AS_reducible\n"
			"\n";
		std::exit(EXIT_SUCCESS);
	} else if (argc == 2) {
		begin = std::stoi(argv[1]);
		end = std::stoi(argv[1]);
	} else if (argc == 3) {
		begin = std::stoi(argv[1]);
		end = std::stoi(argv[2]);
	} else {
		begin = 2;
		end = 10;
	}

	if (begin >= 32 || end >= 32 || begin <= 1) {
		std::cout << " Womp womp try again\n";
		std::exit(EXIT_FAILURE);
	}

	std::cout << std::format("{:>4} {:>12} {:>8} {:>11}\n", "base", "total", "length", "time(s)");
	std::cout << std::string(37, '-') << '\n';
	for (std::size_t base = begin; base <= end; base++) {
		std::chrono::time_point start = std::chrono::steady_clock::now();

		const std::vector<std::size_t> counts = count_parallel_backtrack(base);
		const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

		const std::size_t total = std::reduce(counts.begin(), counts.end());

		std::string counts_str = "[";
		for (std::size_t i = 0; i < counts.size(); ++i)
			counts_str += std::format("{}{}", counts[i], i + 1 < counts.size() ? ", " : "");
		counts_str += "]";

		std::cout << std::format("{:>4} {:>12} {:>8} {:>11.3f}  {}\n", base, total, counts.size() - 1, elapsed, counts_str);
	}

	return 0;
}
