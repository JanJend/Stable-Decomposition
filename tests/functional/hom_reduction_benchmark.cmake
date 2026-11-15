# Run with cmake -DBENCHMARK=/absolute/path/to/hom_reduction_benchmark -P this_file
if(NOT EXISTS "${BENCHMARK}")
    message(FATAL_ERROR "Set BENCHMARK to the built hom_reduction_benchmark executable")
endif()
set(work "${CMAKE_CURRENT_BINARY_DIR}/hom-reduction-fixtures")
file(MAKE_DIRECTORY "${work}")
set(input "${work}/incomparable free generators.scc")
# End_0 is diagonal (dimension 2). At shift 1 all four matrix units exist;
# Q has dimension 2 and one generator per side, while full Hom needs two.
file(WRITE "${input}" "scc2020\n2\n0 2 0\n0 1 ;\n1 0 ;\n")
execute_process(COMMAND "${BENCHMARK}" "${input}" 0 0.5 --adaptive
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE errors TIMEOUT 60)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Benchmark failed: ${output}\n${errors}")
endif()
foreach(expected "End_0_dimension=2" "End_d_basis\t2\t2\t" "End_d_basis\t4\t4\t"
        "quotient_basis\t0\t0\t" "quotient_basis\t2\t2\t"
        "adaptive_phase\tpasses\tcandidate_visits\tretained_visits\treduction_seconds")
    string(FIND "${output}" "${expected}" position)
    if(position EQUAL -1)
        message(FATAL_ERROR "Missing '${expected}': ${output}")
    endif()
endforeach()
foreach(method sparse_nilpotent complement_orbit combined)
    foreach(expected "quotient_${method}\t0\t0\t" "quotient_${method}\t1\t1\t"
            "full_${method}\t1\t1\t" "full_${method}\t2\t2\t")
        string(FIND "${output}" "${expected}" position)
        if(position EQUAL -1)
            message(FATAL_ERROR "Missing '${expected}': ${output}")
        endif()
    endforeach()
endforeach()
# Percentage scaling must use relation degrees too, and the larger box side.
file(WRITE "${input}" "scc2020\n2\n1 2 0\n5 2 ; 0\n0 0 ;\n1 0 ;\n")
execute_process(COMMAND "${BENCHMARK}" "${input}" 10 --percent
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE errors TIMEOUT 30)
if(NOT result EQUAL 0 OR NOT output MATCHES "bounding_box_extent=5 epsilon_units=percent"
        OR NOT output MATCHES "epsilon=0.5 shift=1")
    message(FATAL_ERROR "Percentage epsilon used the wrong bounding box: ${output}\n${errors}")
endif()
# A generator killed by the shift must disappear even on the no-quotient route:
# these sizes count module maps, not distinct lifts of the zero map.
file(WRITE "${input}" "scc2020\n2\n1 1 0\n0.5 0.5 ; 0\n0 0 ;\n")
execute_process(COMMAND "${BENCHMARK}" "${input}" 1
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE errors TIMEOUT 30)
if(NOT result EQUAL 0 OR NOT output MATCHES "End_d_basis\t0\t0\t"
        OR NOT output MATCHES "full_combined\t0\t0\t")
    message(FATAL_ERROR "Zero Hom space was miscounted: ${output}\n${errors}")
endif()
foreach(epsilon -1 nan inf 1e308 0.5junk)
    execute_process(COMMAND "${BENCHMARK}" "${input}" "${epsilon}"
        RESULT_VARIABLE result OUTPUT_QUIET ERROR_QUIET)
    if(NOT result EQUAL 1)
        message(FATAL_ERROR "Invalid epsilon '${epsilon}' was accepted")
    endif()
endforeach()
message(STATUS "Passed Hom dimensions, both generator routes, adaptive run, percentage scaling, zero Hom, and epsilon validation")
