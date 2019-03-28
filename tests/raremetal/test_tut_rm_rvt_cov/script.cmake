# This script is run as a CMake script; therefore the binary directory may not be the same as for the original build
# Use CMake variables to pass in context from outside the script


set(BASE_INPUT_DIR ${BASE_TEST_FOLDER}/inputs)

# Delete old outputs before every single run
set(BASE_OUTPUT_DIR ${BASE_TEST_FOLDER}/output)

if(EXISTS ${BASE_OUTPUT_DIR})
    file(REMOVE_RECURSE ${BASE_OUTPUT_DIR})
    message("Existing output directory removed- ${BASE_OUTPUT_DIR}")
endif()

file(MAKE_DIRECTORY ${BASE_OUTPUT_DIR})


# The input files are relative to the root test directory
execute_process(
        COMMAND ${EXEC_PATH}
            --summaryFiles ${BASE_INPUT_DIR}/rvtests.summaryfiles
            --covFiles ${BASE_INPUT_DIR}/rvtests.covfiles
            --groupFile ${BASE_INPUT_DIR}/group.file
            --burden
            --prefix ${BASE_OUTPUT_DIR}/COMBINED.rvtests.QT1 -
        WORKING_DIRECTORY ${BASE_TEST_FOLDER}
        RESULT_VARIABLE rm_exit_code
        OUTPUT_VARIABLE out_text
        ERROR_VARIABLE out_text)

if(rm_exit_code)
    message("Result exit code: ${rm_exit_code} - ${EXEC_PATH}")
    message(FATAL_ERROR "An error was encountered:\n ${out_text}")
endif()

execute_process(
        COMMAND ${EXEC_PATH}
            --summaryFiles ${BASE_INPUT_DIR}/rmw.summaryfiles
            --covFiles ${BASE_INPUT_DIR}/rmw.covfiles
            --groupFile ${BASE_INPUT_DIR}/group.file
            --burden
            --prefix ${BASE_OUTPUT_DIR}/COMBINED.rmw.QT1 -
        WORKING_DIRECTORY ${BASE_TEST_FOLDER}
        RESULT_VARIABLE rm_exit_code
        OUTPUT_VARIABLE out_text
        ERROR_VARIABLE out_text)

if(rm_exit_code)
    message("Result exit code: ${rm_exit_code} - ${EXEC_PATH}")
    message(FATAL_ERROR "An error was encountered:\n ${out_text}")
endif()

# diff -q COMBINED.rvtests.QT1.meta.burden.results COMBINED.rmw.QT1.meta.burden.results

# Find all output files
# file(GLOB files_list
        # LIST_DIRECTORIES false
        # RELATIVE ${BASE_TEST_FOLDER}/expected
        # ${BASE_TEST_FOLDER}/expected/*.results ${BASE_TEST_FOLDER}/expected/*.tbl)

## Compare the two sets of meta-analysis results to the expected outputs
# foreach(file ${files_list})
set(rvt_b_res COMBINED.rvtests.QT1.meta.burden.results)
set(rmw_b_res COMBINED.rmw.QT1.meta.burden.results)

if(EXISTS ${BASE_OUTPUT_DIR}/${rvt_b_res} AND EXISTS ${BASE_OUTPUT_DIR}/${rmw_b_res})
    message("filename is ${BASE_OUTPUT_DIR}/${file}")
    execute_process(COMMAND diff -q ${BASE_OUTPUT_DIR}/${rvt_b_res} ${BASE_OUTPUT_DIR}/${rmw_b_res} RESULT_VARIABLE diff_exit_code)
    if(diff_exit_code)
        message(FATAL_ERROR "Could not replicate expected calculation results for ${BASE_OUTPUT_DIR}/COMBINED.rvtests.QT1.meta.burden.results and ${BASE_OUTPUT_DIR}/COMBINED.rmw.QT1.meta.burden.results.")
    endif()

else()
    if(NOT EXISTS ${BASE_OUTPUT_DIR}/${rvt_b_res})
        message(FATAL_ERROR "The expected output file ${rvt_b_res} was not generated")
    endif()
    if(NOT EXISTS ${BASE_OUTPUT_DIR}/${rmw_b_res})
        message(FATAL_ERROR "The expected output file ${rmw_b_res} was not generated")
    endif()
endif()
# endforeach()
