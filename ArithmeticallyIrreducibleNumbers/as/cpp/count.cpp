#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <cassert>
#include <chrono>
#include <format>
#include <iostream>
#include <numeric>
#include <ranges>
#include <string>
#include <vector>
#include <atomic>
#include <omp.h>

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
std::vector<std::size_t> count_serial(const std::size_t base, bool conjecture_L = false) {
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

// Parallelisation idea:
// partition the tree so no branch is relatively large,
// each branch is counted using the above Back-Tracking approach,
// separate branches can be computed in parallel

// based on finding the tree in 'base' (base-1) / 8 + 1
// then scaling by a factor of 8 to find the sub-tree
// These are the roots of very expensive branches (heuristic based on the base-17 timings)
// split these into more sub-branches, 
// add another split if the last digit is a multiple of 4
std::pair<std::vector<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>>, std::vector<std::size_t>> partition_tree(
    const size_t base, 
    const bool conjecture_L = false) {
  std::vector<std::size_t> counts(find_upper_bound_L(base, conjecture_L)+2, 0);
  std::vector<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>> branches;

  // organise digits that are 
  //  1) multiples of 8
  //  2) mutliples of 4 but not 8
  //  3) neither
  std::vector<std::size_t> digits_1, digits_2, digits_3;
  for (std::size_t i = 1; i < base; i++) {
    if (i % 8) {
      if (i % 4) {
        digits_3.push_back(i);
      } else {
        digits_2.push_back(i);
      }
    } else {
      digits_1.push_back(i);
    }
  }

  // similar to the original back-tracking algorithm
  const auto pseudo_base = (base - 1) / 8 + 1;
  const auto pseudo_L_upper = find_upper_bound_L(pseudo_base, conjecture_L);
  const auto bound = (pseudo_base - 1) * pseudo_L_upper;
  const auto left_bits = (bound + 2 * (pseudo_base - 1)) / 2; // must get to twice top digit
  const auto right_bits = (bound - 1) / 2;
  const auto total_bits = left_bits + 1 + right_bits;

  std::vector<std::size_t> template_prefix;
  template_prefix.reserve(pseudo_L_upper);
  std::vector<std::size_t> single_end_digit(1, 0);

  std::vector<boost::dynamic_bitset<>> list_of_disallowed_masks;
  list_of_disallowed_masks.reserve(pseudo_L_upper+2);
  boost::dynamic_bitset<> mask(total_bits);
  mask.set(right_bits);
  list_of_disallowed_masks.push_back(mask);

  std::vector<std::size_t> stack;
  stack.reserve((pseudo_base - 1) * pseudo_L_upper);
  for (std::size_t i = 1; i < pseudo_base; i++) {
    stack.push_back(i);
  }

  while (!stack.empty()) {
    std::size_t digit = stack.back();
    stack.pop_back();
    if (!digit) { 
      list_of_disallowed_masks.pop_back();
      template_prefix.pop_back();
      continue;
    }
    stack.emplace_back(0);
    template_prefix.push_back(8 * digit);
    counts[template_prefix.size()] += 1;

    auto mask = list_of_disallowed_masks.back();
    mask = (mask << digit) | (mask >> digit);
    mask.set(right_bits);
    list_of_disallowed_masks.push_back(mask);

    for (std::size_t i = 1; i < pseudo_base; i++){
      if (!mask.test(right_bits + i)){
        stack.emplace_back(i);
      }
    }
    
    // 
    for (auto i : digits_2) {
      auto prefix = template_prefix;
      prefix.push_back(i);
      counts[prefix.size()] += 1;

      // all multiples of 8, but separately because expensive
      for (auto j : digits_1) {
        single_end_digit[0] = j;
        branches.emplace_back(prefix, single_end_digit);
      }

      // only valid multiples of 4, but separately because expensive
      for (auto j : digits_2) {
        if (j == i) {
          continue;
        }
        auto a = (i + j) / 8;
        if (mask.test(right_bits + a)) {
          continue;
        }
        std::size_t b = (i < j) ? (j - i) / 8 : (i - j) / 8;  // abs(i-j)/8
        if (mask.test(right_bits + b)) {
          continue;
        }
        single_end_digit[0] = j;
        branches.emplace_back(prefix, single_end_digit);
      }

      // all other digits, together because cheap
      branches.emplace_back(prefix, digits_3);
    }
    
    // all other digits, together because cheap
    branches.emplace_back(template_prefix, digits_3);

  }

  // start with a non multiple of 8, multiple of 4
  std::vector<std::size_t> prefix(1, 0);
  for (auto i : digits_2) {
    prefix[0] = i; 
    counts[1] += 1;

    // all multiples of 8, but separately becuase expensive
    std::vector<std::size_t> single_end_digit(1, 0);
    for (auto j : digits_1) {
      single_end_digit[0] = j;
      branches.emplace_back(prefix, single_end_digit);
    }

    // only valid multiples of 4, but separately because expensive
    for (auto j : digits_2) {
      if (j == i) {
        continue;
      }
      single_end_digit[0] = j;
      branches.emplace_back(prefix, single_end_digit);
    }

    // all other digits, together because cheap
    branches.emplace_back(prefix, digits_3);
  }
  // finally, odd digits together and even digits separately
  prefix.pop_back();
  std::vector<size_t> even_end_digit;
  even_end_digit.push_back(0);
  std::vector<size_t> odd_end_digits; 
  for (auto i : digits_3) {
    if (i%2) {
      odd_end_digits.push_back(i);
      continue;
    }
    even_end_digit[0] = i;
    branches.emplace_back(prefix, even_end_digit);
  }
  branches.emplace_back(prefix, odd_end_digits);

  // count the branches not yet counted
  for (auto [prefix, end_digits] : branches) {
    counts[prefix.size()+1] += end_digits.size();
  }
  return std::pair{branches, counts};
}

//  count each branch based on the serial back-tracking algorithm
//  we count some small related branches together
//  eg. prefix = [a, b, c], end_digits = [d, e, f]
//  encodes the branches [a, b, c, d], [a, b, c, e] and [a, b, c, f]
std::vector<std::size_t> count_branch(
    std::size_t base,
    std::size_t L_upper,
    std::size_t right_bits,
    std::size_t total_bits,
    std::pair<std::vector<std::size_t>, std::vector<std::size_t>> branch) {
  
  std::vector<boost::dynamic_bitset<>> list_of_disallowed_masks;
  list_of_disallowed_masks.reserve(L_upper+2);
  boost::dynamic_bitset mask(total_bits, 0);
  list_of_disallowed_masks.push_back(mask);
  mask.set(right_bits);
  list_of_disallowed_masks.push_back(mask);
  std::vector<std::size_t> counts(L_upper+2, 0);
  std::vector<std::size_t> stack;
  stack.reserve((base-1) * L_upper ); // probably smaller but ahhh well
  
  auto [prefix, end_digits] = branch;

  // catch the bitmasks up to the prefix
  for (auto & digit : prefix) {
    auto mask = list_of_disallowed_masks.back();
    mask = (mask << digit) | (mask >> digit);
    mask.set(right_bits);
    list_of_disallowed_masks.push_back(mask);
  }

  //  for each end digit in consideration, 
  //  count the entries of the branch (prefix + end_point) of the tree
  for (auto& end_digit : end_digits) {
    stack.push_back(end_digit);
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
  }
  return counts;
}

void reduce(std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
  assert(a.size() == b.size());
  for (std::size_t i = 0; i < a.size(); i++) {
    a[i] += b[i];
  }
}

std::vector<std::size_t> count_parallel(const size_t base, const bool conjecture_L = false) {
  if (base <= 8) {
    return count_serial(base, conjecture_L);	// fast enough, so we focus on parallelising 'big' bases
  }

  std::chrono::time_point start_tree = std::chrono::steady_clock::now();
  std::cout << std::format("\nCounting AS-irreducible strings in base {}...\n", base); 
  auto result = partition_tree(base, conjecture_L);
  auto branches = std::get<0>(result);
  auto counts = std::get<1>(result);
  std::cout << std::format("Search tree partitioned into {} branches.\n\n", branches.size()); 

  auto basics = find_basics(base, conjecture_L); 
  auto L_upper = std::get<0>(basics);
  auto right_bits = std::get<1>(basics);
  auto total_bits = std::get<2>(basics);
  std::cout << std::format("{:>3}/{}  {:>7}  {:>7}  {:>10}  {}\n", "num", branches.size(), "Time", "time", "found", "branch");

  // count along multiple branches in parallel!
  std::atomic<std::size_t> completed{0};
  const auto total_branches = branches.size();
  # pragma omp parallel
  {
    std::vector<std::size_t> local_counts(L_upper + 2, 0);

    #pragma omp for
    for (auto branch : branches) { 
      std::chrono::time_point start_branch = std::chrono::steady_clock::now();

      auto new_counts = count_branch(base, L_upper, right_bits, total_bits, branch);
      reduce(local_counts, new_counts);
 
      const std::size_t sum = std::accumulate(new_counts.begin(), new_counts.end(), 0);
      const double time_branch = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_branch).count();
      const double time_tree = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_tree).count();
      
      auto local_completed = ++completed;
      std::string line = std::format("{:>3}/{}  {:>7.1f}  {:>7.1f}  {:10}  {}\n", local_completed, total_branches, time_tree, time_branch, sum, branch);

      #pragma omp critical (print) 
      {
        std::cout << line;
      }
    }
    #pragma omp critical (account) 
    {
      reduce(counts, local_counts); 
    }
  }
 
  const double time_tree = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_tree).count();
  std::cout << std::format("\nSearch complete! {:.3f}\n\n", time_tree);
  while (!counts.back()) {
    counts.pop_back();
  }
  return counts;
}

void printHelp() {
 std::cout <<
R"(Usage:
  prog [options] [bases]


  Options:

    -h --help
      prints this help statement
    
    -cl --conjecture_L
      uses the conjecture L_max=2^(ceil(log_2(base)))-1 to reduce computation


  There are three modes, by including 0, 1 or 2 integer parameters.

  1. One base, computes that base.
  2. Two bases, computes the inclusive range [start, end].
  3. No bases, computes the default range [2, 16].

  A serial calculation is used for 2 <= base <= 16, 
  a parallel calculation is used for 17 <= base <= 32.
  )";
}

int main(int argc, char *argv[]) {
  bool conjecture_L = false;
  std::vector<int> bases;

  for (int i = 1; i < argc; ++i) {
    std::string_view arg(argv[i]);
    if (arg == "-h" || arg == "--help") {
      printHelp();
      return 0;
    }
    else if (arg == "-cl" || arg == "--conjecture_L") {
      conjecture_L = true;
    } 
    else {
      bases.push_back(std::stoi(std::string(arg)));
    }
  }

  std::size_t start = 2, end = 16;
  if (bases.size() == 1) {
    start = end = bases[0];
  }
  else if (bases.size() == 2) {
    start = bases[0];
    end = bases[1];
    if (start > end)
      throw std::runtime_error("start must be <= end");
  } else if (bases.size() > 2) {
    throw std::runtime_error("At most two base arguments allowed");
  }

  if (start >= 32 || end >= 32 || start <= 1) {
    throw std::runtime_error("bases must be within the range [2, 32]\n");
  }

  std::cout << std::format("{:>4} {:>12} {:>8} {:>11}\n", "base", "total", "length", "time(s)");
  std::cout << std::string(38, '-') << '\n';
  for (std::size_t base = start; base <= end; base++) {
    std::chrono::time_point start = std::chrono::steady_clock::now();

    std::vector<std::size_t> counts; 
    if (base <= 16) {
      counts = count_serial(base, conjecture_L);
    } else {
      counts = count_parallel(base, conjecture_L);
    }

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
