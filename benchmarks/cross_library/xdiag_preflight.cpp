// XDiag v0.6.0 fixed-U(1) capacity probe; source-only until build-verified.
// Use the .bonds file produced by export_case.py. Do not present this as a
// benchmark unless exact source version, solver accuracy, input hash,
// affinity, repeats, and resource allocation are independently attested.
#include <xdiag/all.hpp>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cmath>
#include <limits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <sys/resource.h>

using SteadyClock = std::chrono::steady_clock;

static double seconds(SteadyClock::time_point start, SteadyClock::time_point end) {
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
    if (v > std::numeric_limits<uint64_t>::max()) throw std::overflow_error("binomial overflow");
    res = (uint64_t)v;
  }
  return res;
}
int main(int argc, char **argv) try {
  if (argc != 6 && argc != 7) {
    std::cerr << "Usage: xdiag_preflight N m {dimension|basis|solve} max_u1_dim bonds_file [heisenberg-nn|cubic-j1j2-xxz]\\n";
    return 2;
  }
  int64_t n = std::stoll(argv[1]);
  int64_t m = std::stoll(argv[2]);
  std::string mode = argv[3];
  uint64_t cap = std::stoull(argv[4]);
  std::string filename = argv[5];
  std::string model = (argc == 7 ? argv[6] : "heisenberg-nn");
  if (model != "heisenberg-nn" && model != "cubic-j1j2-xxz")
    throw std::invalid_argument("unknown physical model");
  if (n <= 2 || n > 4096 || m < 0 || m > 3 || m > n || cap == 0 ||
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
  auto start = SteadyClock::now();
  Spinhalf block(n, n - m);
  auto tbasis = SteadyClock::now();
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
  struct Bond { int64_t a, b; double jxy, jzz; };
  std::vector<Bond> bonds;
  std::set<std::pair<int64_t,int64_t>> pairs;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream row(line);
    int64_t i, j;
    double jxy=1.0,jzz=1.0;
    if (!(row >> i >> j)) throw std::runtime_error("malformed bond endpoints");
    if (model == "cubic-j1j2-xxz") {
      if (!(row >> jxy >> jzz))
        throw std::runtime_error("missing Jxy/Jzz on XXZ bond");
    }
    std::string excess;
    if (row >> excess) throw std::runtime_error("extra bond columns");
    if (i < 0 || j >= n || i >= j || !std::isfinite(jxy) || !std::isfinite(jzz))
      throw std::runtime_error("invalid bond endpoints/coefficients");
    if (!pairs.emplace(i,j).second) throw std::runtime_error("duplicate bond");
    bonds.push_back({i,j,jxy,jzz});
  }
  if (bonds.empty()) throw std::runtime_error("no bonds");
  OpSum ops;
  for (auto const &b : bonds) {
    if (model == "cubic-j1j2-xxz") {
      // Exchange = (S+_i S-_j + S-_i S+_j)/2 (XDiag convention).
      // SzSz = Sz_i Sz_j. Identical spin-operator normalization to QuSpin
      // with pauli=0; PR16 ED2 uses Jx=Jy=Jxy and Jz=Jzz.
      ops += b.jxy * Op("Exchange", {b.a, b.b});
      ops += b.jzz * Op("SzSz", {b.a, b.b});
    } else {
      ops += "J" * Op("SdotS", {b.a, b.b});
    }
  }
  if (model == "heisenberg-nn") ops["J"] = 1.0;
  auto tsetup = SteadyClock::now();
  double e0 = eigval0(ops, block); // API examples: SdotS, Exchange, SzSz.
  auto tsolve = SteadyClock::now();
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
