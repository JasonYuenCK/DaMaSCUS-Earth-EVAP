if(NOT DEFINED PROGRAM OR NOT DEFINED CONFIG OR NOT DEFINED EXPECTED_TEXT)
  message(FATAL_ERROR "PROGRAM, CONFIG, and EXPECTED_TEXT must be defined")
endif()

execute_process(
  COMMAND "${PROGRAM}" "${CONFIG}"
  RESULT_VARIABLE command_result
  OUTPUT_VARIABLE command_stdout
  ERROR_VARIABLE command_stderr
  TIMEOUT 5)

if("${command_result}" STREQUAL "0")
  message(FATAL_ERROR "Command unexpectedly succeeded")
endif()
if(NOT "${command_result}" MATCHES "^-?[0-9]+$")
  message(FATAL_ERROR "Command did not exit normally: ${command_result}")
endif()

set(command_output "${command_stdout}\n${command_stderr}")
string(FIND "${command_output}" "${EXPECTED_TEXT}" expected_text_position)
if(expected_text_position EQUAL -1)
  message(FATAL_ERROR
    "Command failed for the wrong reason. Expected: ${EXPECTED_TEXT}\n"
    "Actual output:\n${command_output}")
endif()
