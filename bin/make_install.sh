#!/bin/bash -l
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${PROJECT_DIR}/build"

# For single-config generators such as Unix Makefiles/Ninja, optimization is
# selected at configure time via CMAKE_BUILD_TYPE. The --config argument below
# is only meaningful for multi-config generators such as Xcode/Visual Studio,
# so keep both to make the script portable.
cmake -S "${PROJECT_DIR}" -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCODE_COVERAGE=OFF

cmake --build "${BUILD_DIR}" --config Release --parallel "$(sysctl -n hw.ncpu)"
cmake --install "${BUILD_DIR}" --config Release
