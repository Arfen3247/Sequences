#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <format>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <string>
#include <vector>

std::uint32_t find_upper_bound_L(const std::uint32_t base, bool conjecture_L) {
  if (conjecture_L) { // conjectured upper bound 2^n-1, where
                      // n=ceil(log_2(base))
    return (1u << (32 - __builtin_clz(base - 1))) -
           1; // what is this funny function
  }
  return 2 * base - 1; // safe upper bound
}

std::vector<std::size_t> count_serial(const std::size_t base) {
  // digits in this base
  const std::vector<std::uint16_t> digits = [base] {
    std::vector<std::uint16_t> v(base - 1);
    std::iota(v.begin(), v.end(), 1u);
    return v;
  }();

  // bound the maximum length
  std::uint32_t L_upper = find_upper_bound_L(base, false);

  // bound how many bits are needed
  const std::uint16_t T = (base - 1) * L_upper;
  const std::uint16_t lz = (T + base - 1) / 2;
  const std::uint16_t rz = (T - 1) / 2;
  const std::uint16_t num_bits = lz + 1 + rz;
  boost::dynamic_bitset s0(num_bits, 0);
  s0.set(rz);

  // the sums bitmasks, representing the digits not allowed to appear next
  std::vector<boost::dynamic_bitset<>> sums_list;
  sums_list.reserve(L_upper + 1);
  sums_list.push_back(s0);

  // required result
  std::vector<std::uint64_t> counts(L_upper + 2, 0);
  counts.at(1) = base - 1;

  // the stack for the back-tracking algorithm
  std::vector<std::uint16_t> stack = digits;
  stack.reserve((base - 1) * L_upper);
  while (!stack.empty()) {
    std::uint16_t digit = stack.back();
    stack.pop_back();

    // 0 is the signal to backtrack the last digit
    if (!digit) {
      sums_list.pop_back();
      continue;
    }
    stack.emplace_back(0);

    // apply this digit
    boost::dynamic_bitset s = sums_list.back();
    s = (s << digit) | (s >> digit);
    s.set(rz);
    sums_list.push_back(s);

    // extract candidates for next digit
    std::uint64_t found = 0;
    for (std::uint16_t i = 1; i < base; i++) {
      if (!s.test(rz + i)) {
        stack.emplace_back(i);
        ++found;
      }
    }
    counts.at(sums_list.size()) += found;
  }

  // remove trailing zeros
  while (!counts.back()) {
    counts.pop_back();
  }
  return counts;
}

void branch(std::vector<std::uint16_t> &stack,
            std::vector<boost::dynamic_bitset<>> &sums_list,
            std::vector<std::uint64_t> &counts, const std::size_t base,
            const std::uint16_t rz) {
  std::uint16_t digit = stack.back();
  stack.pop_back();

  // 0 is the signal to backtrack the last digit
  if (!digit) {
    sums_list.pop_back();
    return;
  }
  stack.push_back(0);

  // apply this digit
  boost::dynamic_bitset s = sums_list.back();
  s = (s << digit) | (s >> digit);
  s.set(rz);
  sums_list.push_back(s);

  // extract candidates for next digit
  std::uint64_t found = 0;
  for (std::uint16_t i = 1; i < base; i++) {
    if (!s.test(rz + i)) {
      stack.push_back(i);
      ++found;
    }
  }
  counts.at(sums_list.size()) += found;
  // counts[sums_list.size()] += found;
}

void extend_stack(std::vector<std::uint16_t> &omp_out,
                  const std::vector<std::uint16_t> &omp_in) {
  omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end());
}
void reduce_count(std::vector<std::uint64_t> &omp_out,
                  const std::vector<std::uint64_t> &omp_in) {
  for (std::size_t i = 0; i < omp_out.size(); i++) {
    omp_out.at(i) += omp_in.at(i);
  }
}
std::vector<std::size_t> count_parallel(const std::size_t base) {
  // digits in this base
  const std::vector<std::uint16_t> digits = [base] {
    std::vector<std::uint16_t> v(base - 1);
    std::iota(v.begin(), v.end(), 1u);
    return v;
  }();

  // bound the maximum length
  const std::uint32_t L_upper = find_upper_bound_L(base, false);

  // bound how many bits are needed
  const std::uint16_t T = (base - 1) * L_upper;
  const std::uint16_t lz = (T + base - 1) / 2;
  const std::uint16_t rz = (T - 1) / 2;
  const std::uint16_t num_bits = lz + 1 + rz;
  const boost::dynamic_bitset s0 = [num_bits, rz]() {
    boost::dynamic_bitset s0(num_bits, 0);
    s0.set(rz);
    return s0;
  }();

  // the sums bitmasks, representing the digits not allowed to appear next
  std::vector<boost::dynamic_bitset<>> sums_list;
  sums_list.reserve(L_upper + 1);
  sums_list.push_back(s0);
// required result
  std::vector<std::uint64_t> counts(L_upper + 2, 0);
  counts.at(1) = base - 1;

  // the stack for the back-tracking algorithm
  std::vector<std::uint16_t> stack = digits;
  stack.reserve((base - 1) * L_upper);

  // #pragma omp parallel for  reduction(extend :
  // counts) schedule(dynamic)

#pragma omp parallel for firstprivate(sums_list) schedule(dynamic)
  for (std::size_t i = 0; i < stack.size(); i++) {
    std::vector<std::uint16_t> local_stack = {stack.at(i)};
    std::vector<std::uint64_t> local_counts(counts.size(), 0);
    while (!local_stack.empty()) {
      branch(local_stack, sums_list, local_counts, base, rz);
    }

#pragma omp critical
    reduce_count(counts, local_counts);
  }

  //  while (!stack.empty()) {
  //    branch(stack, sums_list, rz, counts, base);
  //  }

  // remove trailing zeros
  while (!counts.back()) {
    counts.pop_back();
  }
  return counts;
}

enum mode { serial = 0, parallel = 1 };

std::vector<std::size_t> count(const std::size_t base,
                               const mode mode = serial) {
  switch (mode) { // completely overkill here using switch case but ahh well
  case serial:
    return count_serial(base);
  case parallel:
    return count_parallel(base);
  default:
    return count_serial(base);
  }
}

int main(int argc, char *argv[]) {

  std::optional<std::size_t> begin, end;
  std::optional<mode> exec_mode;

  const std::vector<std::string> args = [argc, argv]() {
    std::vector<std::string> args;
    for (std::size_t i = 1; i < static_cast<std::size_t>(argc); i++) {
      args.emplace_back(argv[i]);
    }
    return args;
  }();

  const auto info = []() {
    std::cout << "The following shows examples of how this binary is designed "
                 "to be used\n"
                 "\n"
                 "Examples on how to run this binary\n"
                 "To count only in base 5\n"
                 "    ./AS_reducible -begin 5\n"
                 "    ./AS_reducible -b 5\n"
                 "To count in range from 5 to 10\n"
                 "    ./AS_reducible -begin 5 -end 10\n"
                 "    ./AS_reducible -b 5 -e 10\n"
                 "To select mode either serial or parallel\n"
                 "    ./AS_reducible -mode parallel -b 10\n"
                 "    ./AS_reducible -mode serial -b 10\n"
                 "    ./AS_reducible -m p -b 10\n"
                 "    ./AS_reducible -m s -b 10\n"
                 "\n";
  };

  // for (auto it =  args.begin(); it != args.end(); it++ ) {
  for (std::size_t i = 0; i < args.size(); i++) {
    std::string arg = args.at(i);
    std::cout << arg << std::endl;

    if (arg == "-b" || arg == "-begin") {
      begin = std::stoi(args.at(i + 1));
      i += 1;
      continue;
    }

    if (arg == "-e" || arg == "-end") {
      end = std::stoi(args.at(i + 1));
      i += 1;
      continue;
    }

    if (arg == "-m" || arg == "-mode") {
      if (args.at(i + 1) == "parallel" || args.at(i + 1) == "p") {
        exec_mode = parallel;
      } else if (args.at(i + 1) == "serial" || args.at(i + 1) == "s") {
        exec_mode = serial;
      }
      i += 1;
      continue;
    }

    if ((std::string(arg) == "help" || std::string(arg) == "-help" ||
         std::string(arg) == "-h") ||
        std::string(arg) == "h") {
      info();
      std::exit(EXIT_SUCCESS);
    }
  }

  if (!begin && !end && !exec_mode) {
    info();
    std::exit(EXIT_FAILURE);
  }

  if (begin && !end) {
    end = begin.value();
  }
  if (!exec_mode) {
    exec_mode = serial;
  }

  std::cout << std::format("{:>4} {:>12} {:>8} {:>11}\n", "base", "total",
                           "length", "time(s)");
  std::cout << std::string(37, '-') << '\n';
  for (std::size_t base = begin.value(); base <= end; base++) {
    std::chrono::time_point start = std::chrono::steady_clock::now();

    const std::vector<std::size_t> counts = count(base, exec_mode.value());
    const double elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
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
