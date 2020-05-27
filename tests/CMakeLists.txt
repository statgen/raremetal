# Currently only raremetalworker is tested here.
# raremetal has a new Catch2 test suite. We'll port these remaining tests to raremetalworker down the road,
# assuming we don't completely deprecate raremetalworker.

set(RMW_BIN ${CMAKE_BINARY_DIR}/raremetalworker)


##### Definitions and helpers
function(copy_test_scaffold
        BASE_TEST_DESTINATION
        INPUT_SRC_DIR # where to find the inputs that will be copied over. May leave blank
        EXPECTED_SRC_DIR)  # where to find the expected outputs; may leave blank
    ## Prepare the basic directory structure (in the build folder) for tests to run
    # Each test lives under a directory with the same name as the test, and has a base folder with
    #   `inputs`, `expected`, and `outputs`. (the first two are copied over into `build/` as needed, the latter is
    #   wiped and regenerated every time the test runs)
    file(MAKE_DIRECTORY ${BASE_TEST_DESTINATION})

    if(INPUT_SRC_DIR)
        file(MAKE_DIRECTORY ${BASE_TEST_DESTINATION}/inputs)
        file(COPY ${INPUT_SRC_DIR}
                DESTINATION ${BASE_TEST_DESTINATION}/inputs)
    endif()


    if(EXPECTED_SRC_DIR)
        file(MAKE_DIRECTORY ${BASE_TEST_DESTINATION}/expected)
        file(COPY ${EXPECTED_SRC_DIR}
                DESTINATION ${BASE_TEST_DESTINATION}/expected)
    endif()



endfunction()

function(run_boilerplate
         BASE_TEST_DESTINATION
         TEST_NAME
         SCRIPT_PATH
         EXEC_PATH)
    # SCRIPT_PATH should be relative to the top-level tests folder
    # Define a single unit test of the specified name, telling it where to find the built test binary
    add_test(NAME ${TEST_NAME}
             COMMAND ${CMAKE_COMMAND}
             -DEXEC_PATH=${EXEC_PATH}
             -DBASE_TEST_FOLDER=${BASE_TEST_DESTINATION}
             -P ${CMAKE_CURRENT_SOURCE_DIR}/${SCRIPT_PATH})
endfunction()


###########################################################################
#  Tests based on the RAREMETAL tutorial dataset as of this repo version  #
###########################################################################

########## RAREMETALWORKER

# When kinship options are specified, RMW uses that information as part of the covariance calcs
set(BASE_TEST_FOLDER ${CMAKE_CURRENT_BINARY_DIR}/test_tut_rmw_kin)
copy_test_scaffold(
        ${BASE_TEST_FOLDER}
        datasets/raremetal_tutorial/
        raremetalworker/test_tut_rmw_kin/expected/)
run_boilerplate(
        ${BASE_TEST_FOLDER}
        test_tut_rmw_kin
        raremetalworker/test_tut_rmw_kin/script.cmake
        ${RMW_BIN})

# When no kinship options are specified, covariance is calculated according to a single fixed sigma term
set(BASE_TEST_FOLDER ${CMAKE_CURRENT_BINARY_DIR}/test_tut_rmw_nokin)
copy_test_scaffold(
        ${BASE_TEST_FOLDER}
        datasets/raremetal_tutorial/
        raremetalworker/test_tut_rmw_nokin/expected/)
run_boilerplate(
        ${BASE_TEST_FOLDER}
        test_tut_rmw_nokin
        raremetalworker/test_tut_rmw_nokin/script.cmake
        ${RMW_BIN})
