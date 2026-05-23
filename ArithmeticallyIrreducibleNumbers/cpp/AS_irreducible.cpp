#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <cassert>
#include <chrono>
#include <format>
#include <iostream>
#include <numeric>
#include <ranges>
#include <string>
#include <vector>

std::vector<std::size_t> count_AS_irreducible(const std::size_t base) {

  const std::vector<std::size_t> digits = [base]() {
    std::vector<std::size_t> digits(base - 1);
    std::iota(digits.begin(), digits.end(), 1);
    return digits;
  }();
  const std::vector<std::size_t> bounds = [base]() {
    std::vector<std::size_t> bounds(2 * base);
    for (std::size_t length_ = 0; length_ < bounds.size(); length_++) {
      bounds.at(length_) = base * (length_ + 1) / 2;
    }
    return bounds;
  }();

  std::vector<std::size_t> counts(2 * base, 0);
  counts.at(1) = base - 1;

  std::vector<std::size_t> arr;

  const std::size_t bitset_length = 2 * base * base;

  // stack will have at most b*L_max states
  typedef std::pair<std::size_t, bool> stackpair;
  std::vector<stackpair> stack(digits.size(), stackpair(0, false));
  for (std::size_t i = 0; i < digits.size(); i++) {
    stack.at(i).first = digits.at(i);
  }

  while (!stack.empty()) {
    std::size_t digit = stack.back().first;
    bool second_pass = stack.back().second;
    stack.pop_back();

    if (second_pass) { // backtrack changes
      arr.pop_back();
      continue;
    }

    stack.emplace_back(digit, true);

    arr.push_back(digit); // enact changes

    // find the digits x such that no suffix of arr+[x] has the AS property,
    // ie. exists signed summation sigma such that sigma(suffix)=0.
    // O(bL) space, O(bL^2) time using bit manipulations (need 2*b*L bits, maybe
    // b*L).

    // NEEDS UPDATING FOR WHEN BASE >= 64
    boost::dynamic_bitset cm(bitset_length, (1 << base) - 2);
    std::size_t bound = bounds.at(arr.size() + 1);
    boost::dynamic_bitset sums(bitset_length, 1);
    sums <<= bound;
    for (std::size_t &y : arr | std::views::reverse) {
      sums = (sums << y) | (sums >> y);
      cm &= ~(sums >> bound);
      if (cm.none()) {
        break;
      }
    }
    if (cm.none()) {
      continue;
    }

    // read off the allowed candidates x in O(b^2) time, add to stack
    for (const std::size_t &x : digits) {
      // if bit that represents x in cm then append to
      if (cm.test(x)) {
        stack.emplace_back(x, false);
        counts.at(arr.size() + 1) += 1;
      }
    }
  }

  while (!counts.back()) {
    counts.pop_back();
  }

  return counts;
}

int main(int argc, char *argv[]) {

  std::size_t begin, end;

  if (argc == 2) {
    begin = std::stoi(argv[1]);
    end = std::stoi(argv[1]);
  } else if (argc == 3) {
    begin = std::stoi(argv[1]);
    end = std::stoi(argv[2]);
  } else if (argc == 3) {
    begin = std::stoi(argv[1]);
    end = std::stoi(argv[2]);
  } else {
    begin = 2;
    end = 11;
  }

  // below formatting and timing is the work of an llm, is due some reviewing but on first glance it looks good
  std::cout << std::format("{:>4} {:>12} {:>8} {:>11}\n", "base", "total",
                           "length", "time(s)");
  std::cout << std::string(37, '-') << '\n';

  for (std::size_t base = begin; base <= end; base++) {
    auto st = std::chrono::steady_clock::now();

    const std::vector<std::size_t> counts = count_AS_irreducible(base);
    double elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - st)
            .count();

    const std::size_t total = std::reduce(counts.begin(), counts.end());

    std::string counts_str = "[";
    for (std::size_t i = 0; i < counts.size(); ++i)
      counts_str +=
          std::format("{}{}", counts[i], i + 1 < counts.size() ? ", " : "");
    counts_str += "]";

    std::cout << std::format("{:>4} {:>12} {:>8} {:>11.3f}  {}\n", base, total,
                             counts.size() - 1, elapsed, counts_str);
  }

  return 0;
}
