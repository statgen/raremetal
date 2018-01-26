# This script is run as a CMake script; therefore the binary directory may not be the same as for the original build
# Use CMake variables to pass in context from outside the script


# Delete old outputs before every single run
set(BASE_OUTPUT_DIR ${BASE_TEST_FOLDER}/output)

if(EXISTS ${BASE_OUTPUT_DIR})
    file(REMOVE_RECURSE ${BASE_OUTPUT_DIR})
    message("Existing output directory removed- ${BASE_OUTPUT_DIR}")
endif()

file(MAKE_DIRECTORY ${BASE_OUTPUT_DIR})

execute_process(
        COMMAND ${RMW_BIN}
            --ped ${CMAKE_BINARY_DIR}/raremetal_tutorial/example1.ped
            --dat ${CMAKE_BINARY_DIR}/raremetal_tutorial/example1.dat
            --vcf ${CMAKE_BINARY_DIR}/raremetal_tutorial/example1.vcf.gz
            --traitName QT1  --inverseNormal --makeResiduals --kinSave --kinGeno --prefix STUDY1
        WORKING_DIRECTORY ${BASE_OUTPUT_DIR}
        RESULT_VARIABLE rmw_exit_code
        OUTPUT_VARIABLE out_text
        ERROR_VARIABLE out_text)

if(rmw_exit_code)
    message("Result exit code: ${rmw_exit_code} - ${RMW_BIN}")
    message(FATAL_ERROR "An error was encountered:\n ${out_text}")
endif()

execute_process(
        COMMAND ${RMW_BIN}
        --ped ${CMAKE_CURRENT_BINARY_DIR}/raremetal_tutorial/example2.ped
        --dat ${CMAKE_CURRENT_BINARY_DIR}/raremetal_tutorial/example2.dat
        --vcf ${CMAKE_CURRENT_BINARY_DIR}/raremetal_tutorial/example2.vcf.gz
        --traitName QT1  --inverseNormal --makeResiduals --kinSave --kinGeno --prefix STUDY2
        WORKING_DIRECTORY ${BASE_OUTPUT_DIR}
        RESULT_VARIABLE rmw_exit_code
        OUTPUT_VARIABLE out_text
        ERROR_VARIABLE out_text)

if(rmw_exit_code)
    message("Result exit code: ${rmw_exit_code} - ${RMW_BIN}")
    message(FATAL_ERROR "an error was encountered: ${out_text}")
endif()

# Find all output files
file(GLOB files_list
        LIST_DIRECTORIES false
        RELATIVE ${BASE_TEST_FOLDER}/expected
        ${BASE_TEST_FOLDER}/expected/*.txt)

## Compare the two sets of meta-analysis results to the expected outputs
foreach(file ${files_list})
    if(EXISTS ${BASE_OUTPUT_DIR}/${file})
        message("filename is ${BASE_OUTPUT_DIR}/${file}")
        execute_process(COMMAND diff -q ${BASE_TEST_FOLDER}/expected/${file} ${BASE_OUTPUT_DIR}/${file} RESULT_VARIABLE diff_exit_code)
        if(diff_exit_code)
            message(FATAL_ERROR "Could not replicate expected calculation results for ${file}.")
        endif()

    else()
        message(FATAL_ERROR "The expected output file ${file} was not generated")
    endif()
endforeach()
