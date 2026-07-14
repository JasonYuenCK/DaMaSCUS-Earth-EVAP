if(NOT DEFINED OBSCURA_SOURCE_DIR)
  message(FATAL_ERROR "OBSCURA_SOURCE_DIR is required")
endif()

set(CRYSTAL_SOURCE "${OBSCURA_SOURCE_DIR}/src/Direct_Detection_Crystal.cpp")
if(NOT EXISTS "${CRYSTAL_SOURCE}")
  message(FATAL_ERROR "obscura crystal source not found: ${CRYSTAL_SOURCE}")
endif()

file(READ "${CRYSTAL_SOURCE}" CONTENTS)
if(CONTENTS MATCHES "using_Q_threshold\\(false\\).+using_Q_bins\\(false\\)")
  return()
endif()

set(PATCHED "${CONTENTS}")
string(REPLACE
  ": DM_Detector(\"Crystal experiment\", gram * year, \"Electrons\"), target_crystal(Crystal(\"Si\"))"
  ": DM_Detector(\"Crystal experiment\", gram * year, \"Electrons\"), target_crystal(Crystal(\"Si\")), Q_threshold(0), using_Q_threshold(false), using_Q_bins(false)"
  PATCHED "${PATCHED}")
string(REPLACE
  ": DM_Detector(label, expo, \"Electrons\"), target_crystal(Crystal(crys))"
  ": DM_Detector(label, expo, \"Electrons\"), target_crystal(Crystal(crys)), Q_threshold(0), using_Q_threshold(false), using_Q_bins(false)"
  PATCHED "${PATCHED}")

if(PATCHED STREQUAL CONTENTS)
  message(FATAL_ERROR
    "unsupported obscura Direct_Detection_Crystal.cpp: constructor pattern changed")
endif()

file(WRITE "${CRYSTAL_SOURCE}" "${PATCHED}")
message(STATUS "Patched obscura crystal detector flag initialization")
