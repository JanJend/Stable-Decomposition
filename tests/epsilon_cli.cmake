function(check_success expected)
    execute_process(COMMAND "${PRUNING}" "${INPUT}" --no-output ${ARGN}
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 30)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "Pruning failed (${result}): ${output}\n${error}")
    endif()
    string(FIND "${output}" "${expected}" position)
    if(position EQUAL -1)
        message(FATAL_ERROR "Missing expected output '${expected}': ${output}")
    endif()
endfunction()

# test1.scc has degree extent 4, hence epsilon=0.04 and shift=0.08.
check_success("(epsilon=0.04, shift=0.08)")
check_success("(epsilon=0.25, shift=0.5)" --epsilon 0.25)
check_success("(epsilon=0.25, shift=0.5)" --delta 0.25)
check_success("(epsilon=0, shift=0)" --epsilon 0)
foreach(value IN ITEMS -1 nan inf 1e308 invalid)
    execute_process(COMMAND "${PRUNING}" "${INPUT}" --epsilon "${value}" --no-output
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 30)
    if(NOT result EQUAL 1)
        message(FATAL_ERROR "Invalid epsilon '${value}' was not rejected cleanly: ${result}")
    endif()
endforeach()
execute_process(COMMAND "${PRUNING}" "${INPUT}" --epsilon
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 30)
if(NOT result EQUAL 1)
    message(FATAL_ERROR "Missing epsilon value was not rejected cleanly: ${result}")
endif()
