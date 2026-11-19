// XDiag v0.6.0 fixed-U(1) capacity probe; source-only until build-verified.
// Use the .bonds file produced by export_case.py. Do not present this as a
// benchmark unless exact source version, solver accuracy, input hash,
// affinity, repeats, and resource allocation are independently attested.
#include <xdiag/all.hpp>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <sys/resource.h>

using clock_t = std::chrono::steady_clock;

static double seconds(clock_t::time_point start, clock_t::time_point end) {
  return std::chrono::duration<double>(end - start).count();
}
static long peak_kib() {
  struct rusage ru {};
  if (getrusage(RUSAGE_SELF, &ru)) return -1;
  return ru.ru_maxrss;  // KiB on Linux only
}
static uint64_t choose_small(int64_t n, int64_t k) {
  if (k < 0 || k > n) throw std::invalid_argument("invalid ndown");
  k = std::min(k, n - k);
  uint64_t res = 1;
  for (int64_t i = 1; i <= k; ++i) {
    // For this pilot k is 2 or 3; keep arithmetic exact and fail closed.
    __uint128_t v = (__uint128_t)res * (n - k + i);
    v /= i;
    if (v > UINT64_MAX) throw std::overflow_error("binomial overflow");
    res = (uint64_t)v;
  }
  return res;
}
int main(int argc, char **argv) try {
  if (argc != 6) {
    std::cerr << "Usage: xdiag_preflight N m {dimension|basis|solve} max_u1_dim bonds_file\n";
    return 2;
  }
  int64_t n = std::stoll(argv[1]);
  int64_t m = std::stoll(argv[2]);
  std::string mode = argv[3];
  uint64_t cap = std::stoull(argv[4]);
  std::string filename = argv[5];
  if (n <= 2 || m < 0 || m > 3 || m > n || cap == 0 ||
      (mode != "dimension" && mode != "basis" && mode != "solve")) {
    throw std::invalid_argument("unsupported N/m/mode/cap");
  }
  const auto dim = choose_small(n, m);
  if (mode == "dimension") {
    std::cout << "{\"status\":\"DIMENSIONS_ONLY\",\"N\":" << n
              << ",\"m\":" << m << ",\"u1_dimension\":" << dim << "}\n";
    return 0;
  }
  if (dim > cap) {
    std::cout << "{\"status\":\"SKIPPED_RESOURCE_CAP\",\"u1_dimension\":"
              << dim << ",\"max_u1_dim\":" << cap << "}\n";
    return 0;
  }

  // This is intentionally only the U(1) block: NO false Gamma/Oh claim.
  using namespace xdiag;
  auto start = clock_t::now();
  Spinhalf block(n, n - m);
  auto tbasis = clock_t::now();
  if ((uint64_t)block.size() != dim) {
    throw std::runtime_error("XDiag U1 dimension mismatch");
  }
  if (mode == "basis") {
    std::cout << "{\"status\":\"BASIS_ONLY\",\"N\":" << n
              << ",\"m\":" << m << ",\"dimension\":" << block.size()
              << ",\"basis_seconds\":" << seconds(start, tbasis)
              << ",\"peak_rss_kib\":" << peak_kib() << "}\n";
    return 0;
  }

  // Single source of physical bonds: canonical export_case.py bond file.
  std::ifstream in(filename);
  if (!in) throw std::runtime_error("cannot open bond file");
  std::set<std::pair<int64_t,int64_t>> bonds;
  int64_t i, j;
  while (in >> i >> j) {
    if (i < 0 || j >= n || i >= j) {
      throw std::runtime_error("invalid bond or order");
    }
    if (!bonds.emplace(i,j).second) {
      throw std::runtime_error("duplicate bond");
    }
  }
  if (!in.eof()) throw std::runtime_error("malformed bond file");
  if (bonds.empty()) throw std::runtime_error("no bonds");
  OpSum ops;
  for (auto [a, b] : bonds) ops += "J" * Op("SdotS", {a, b});
  ops["J"] = 1.0;
  auto tsetup = clock_t::now();
  double e0 = eigval0(ops, block); // API verified against XDiag example.
  auto tsolve = clock_t::now();
  std::cout << "{\"status\":\"ENERGY_ONLY_RESIDUAL_UNCHECKED\",\"N\":" << n
            << ",\"m\":" << m << ",\"dimension\":" << block.size()
            << ",\"bonds\":" << bonds.size() << ",\"E0\":" << e0
            << ",\"basis_seconds\":" << seconds(start, tbasis)
            << ",\"operator_spec_seconds\":" << seconds(tbasis,tsetup)
            << ",\"solve_seconds\":" << seconds(tsetup,tsolve)
            << ",\"time_to_eigenvalue_seconds\":" << seconds(start,tsolve)
            << ",\"peak_rss_kib\":" << peak_kib() << "}\n";
  // E0 alone cannot satisfy the publication-grade direct residual gate.
  return 0;
} catch (std::exception const &e) {
  std::cerr << "XDiag pilot error: " << e.what() << "\n";
  return 3;
}
