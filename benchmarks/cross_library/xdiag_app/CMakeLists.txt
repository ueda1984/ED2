cmake_minimum_required(VERSION 3.19)
project(ed2_xdiag_comparator LANGUAGES CXX)

# Compile against a pinned, externally installed XDiag, not an ED2 fork.
# Usage: cmake -S benchmarks/cross_library/xdiag_app -B build-app \
#          -DCMAKE_PREFIX_PATH=/path/to/xdiag/install
find_package(xdiag REQUIRED CONFIG)
add_executable(xdiag_preflight ../xdiag_preflight.cpp)
target_compile_features(xdiag_preflight PRIVATE cxx_std_17)
target_link_libraries(xdiag_preflight PRIVATE xdiag::xdiag)
