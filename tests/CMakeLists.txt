# 1. Test executable
add_test(NAME executable_test
	COMMAND DaMaSCUS-SUN ${PROJECT_SOURCE_DIR}/tests/config_executable_smoke.cfg
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_test(NAME test_Invalid_Snapshot_Interval_Config
	COMMAND ${CMAKE_COMMAND}
		-DPROGRAM=$<TARGET_FILE:DaMaSCUS-SUN>
		-DCONFIG=${PROJECT_SOURCE_DIR}/tests/config_snapshot_invalid_fractional.cfg
		"-DEXPECTED_TEXT='snapshot_interval' must be a positive integer number of seconds."
		-P ${PROJECT_SOURCE_DIR}/cmake/ExpectCommandFailure.cmake
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(test_Invalid_Snapshot_Interval_Config PROPERTIES
	TIMEOUT 8
	LABELS "snapshot;config")

# 2. Unit tests: Dark_Photon.xpp	
add_executable(test_Dark_Photon
	test_Dark_Photon.cpp)
target_include_directories(test_Dark_Photon PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Dark_Photon
	PRIVATE
		lib_damascus_sun
		gtest_main)	#contains the main function
target_compile_options(test_Dark_Photon PUBLIC -Wall -pedantic)
install(TARGETS test_Dark_Photon DESTINATION ${TESTS_DIR})
add_test(NAME test_Dark_Photon COMMAND test_Dark_Photon
	WORKING_DIRECTORY ${TESTS_DIR})

# 3. Unit tests: Data_Generation.xpp	
add_executable(test_Data_Generation
	test_Data_Generation.cpp)
target_include_directories(test_Data_Generation PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Data_Generation
	PRIVATE
		lib_damascus_sun
		gtest)
target_compile_options(test_Data_Generation PUBLIC -Wall -pedantic)
install(TARGETS test_Data_Generation DESTINATION ${TESTS_DIR})
add_test(NAME test_Data_Generation COMMAND test_Data_Generation
	WORKING_DIRECTORY ${TESTS_DIR})

# 4. Unit tests: Parameter_Scan.xpp	
add_executable(test_Parameter_Scan
	test_Parameter_Scan.cpp)
target_include_directories(test_Parameter_Scan PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Parameter_Scan
	PRIVATE
		lib_damascus_sun
		gtest_main)	#contains the main function
target_compile_options(test_Parameter_Scan PUBLIC -Wall -pedantic)
install(TARGETS test_Parameter_Scan DESTINATION ${TESTS_DIR})
add_test(NAME test_Parameter_Scan COMMAND test_Parameter_Scan
	WORKING_DIRECTORY ${TESTS_DIR})

# 5. Unit tests: Reflection_Spectrum.xpp	
add_executable(test_Reflection_Spectrum
	test_Reflection_Spectrum.cpp)
target_include_directories(test_Reflection_Spectrum PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Reflection_Spectrum
	PRIVATE
		lib_damascus_sun
		gtest_main)	#contains the main function
target_compile_options(test_Reflection_Spectrum PUBLIC -Wall -pedantic)
install(TARGETS test_Reflection_Spectrum DESTINATION ${TESTS_DIR})
add_test(NAME test_Reflection_Spectrum COMMAND test_Reflection_Spectrum
	WORKING_DIRECTORY ${TESTS_DIR})

# 6. Unit tests: Simulation_Trajectory.xpp	
add_executable(test_Simulation_Trajectory
	test_Simulation_Trajectory.cpp)
target_include_directories(test_Simulation_Trajectory PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Simulation_Trajectory
	PRIVATE
		lib_damascus_sun
		gtest_main)	#contains the main function
target_compile_options(test_Simulation_Trajectory PUBLIC -Wall -pedantic)
install(TARGETS test_Simulation_Trajectory DESTINATION ${TESTS_DIR})
add_test(NAME test_Simulation_Trajectory COMMAND test_Simulation_Trajectory
	WORKING_DIRECTORY ${TESTS_DIR})

# 7. Unit tests: Simulation_Utilties.xpp	
add_executable(test_Simulation_Utilities
	test_Simulation_Utilities.cpp)
target_include_directories(test_Simulation_Utilities PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Simulation_Utilities
	PRIVATE
		lib_damascus_sun
		gtest_main)	#contains the main function
target_compile_options(test_Simulation_Utilities PUBLIC -Wall -pedantic)
install(TARGETS test_Simulation_Utilities DESTINATION ${TESTS_DIR})
add_test(NAME test_Simulation_Utilities COMMAND test_Simulation_Utilities
	WORKING_DIRECTORY ${TESTS_DIR})

# 8. Unit tests: Solar_Model.xpp	
add_executable(test_Solar_Model
	test_Solar_Model.cpp)
target_include_directories(test_Solar_Model PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Solar_Model
	PRIVATE
		lib_damascus_sun
		gtest_main)	#contains the main function
target_compile_options(test_Solar_Model PUBLIC -Wall -pedantic)
install(TARGETS test_Solar_Model DESTINATION ${TESTS_DIR})
add_test(NAME test_Solar_Model COMMAND test_Solar_Model
	WORKING_DIRECTORY ${TESTS_DIR})

# 9. Snapshot serialization and merge tests (no solar model/runtime data)
add_executable(test_Snapshot_IO
	test_Snapshot_IO.cpp)
target_include_directories(test_Snapshot_IO PRIVATE ${GENERATED_DIR})
target_link_libraries(test_Snapshot_IO
	PRIVATE
		lib_damascus_sun
		gtest_main)
target_compile_options(test_Snapshot_IO PUBLIC -Wall -pedantic)
install(TARGETS test_Snapshot_IO DESTINATION ${TESTS_DIR})
add_test(NAME test_Snapshot_IO COMMAND test_Snapshot_IO
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(test_Snapshot_IO PROPERTIES
	TIMEOUT 10
	LABELS "snapshot;unit")

# 10. Real four-rank heartbeat test with barrier-synchronized local epochs.
add_executable(test_Snapshot_Heartbeat_MPI
	test_Snapshot_Heartbeat_MPI.cpp)
target_include_directories(test_Snapshot_Heartbeat_MPI PRIVATE ${GENERATED_DIR})
target_link_libraries(test_Snapshot_Heartbeat_MPI
	PRIVATE
		lib_damascus_sun
		${MPI_CXX_LIBRARIES})
target_compile_options(test_Snapshot_Heartbeat_MPI PUBLIC -Wall -pedantic)
install(TARGETS test_Snapshot_Heartbeat_MPI DESTINATION ${TESTS_DIR})
if(MPIEXEC_EXECUTABLE
	AND MPIEXEC_NUMPROC_FLAG
	AND DEFINED MPIEXEC_MAX_NUMPROCS
	AND MPIEXEC_MAX_NUMPROCS GREATER_EQUAL 4)
	add_test(NAME test_Snapshot_Heartbeat_MPI_4Ranks
		COMMAND
			${MPIEXEC_EXECUTABLE}
			${MPIEXEC_NUMPROC_FLAG} 4
			${MPIEXEC_PREFLAGS}
			$<TARGET_FILE:test_Snapshot_Heartbeat_MPI>
			${MPIEXEC_POSTFLAGS}
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
	set_tests_properties(test_Snapshot_Heartbeat_MPI_4Ranks PROPERTIES
		PROCESSORS 4
		RUN_SERIAL TRUE
		TIMEOUT 35
		LABELS "snapshot;mpi;slow")
else()
	message(STATUS "Skipping four-rank snapshot heartbeat test: MPI launcher or explicit four-rank capacity is unavailable (MPIEXEC_EXECUTABLE='${MPIEXEC_EXECUTABLE}', MPIEXEC_MAX_NUMPROCS='${MPIEXEC_MAX_NUMPROCS}')")
endif()
