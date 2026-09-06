cmake_minimum_required(VERSION 3.18)
project(AtmosModel LANGUAGES Fortran)

# Option to control OpenACC vs OpenMP build (like USE_OPENACC in the Makefile)
option(USE_OPENACC "Enable OpenACC (NVHPC) instead of OpenMP" OFF)

# ------------------------------------------------------------
# NetCDF (Fortran + C) flags via nf-config / nc-config
# ------------------------------------------------------------
execute_process(
    COMMAND nf-config --fflags
    OUTPUT_VARIABLE NF_FFLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND nf-config --flibs
    OUTPUT_VARIABLE NF_FLIBS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND nc-config --libdir
    OUTPUT_VARIABLE NC_LIBDIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Convert the strings to CMake lists
separate_arguments(NF_FFLAGS)
separate_arguments(NF_FLIBS)

# ------------------------------------------------------------
# MPI
# ------------------------------------------------------------
find_package(MPI REQUIRED Fortran)

# ------------------------------------------------------------
# Source files
# ------------------------------------------------------------
set(MODEL_SOURCES
    module_parameters.f90
    module_types.F90
    module_output.F90
    module_physics.f90
    parallel_timer.f90
    model.F90
)

add_executable(model ${MODEL_SOURCES})

# ------------------------------------------------------------
# Compile options (mirroring the Makefile)
# ------------------------------------------------------------
# Base NetCDF flags
target_compile_options(model PRIVATE ${NF_FFLAGS})

# Choose between OpenACC and OpenMP-style build
if(USE_OPENACC)
    message(STATUS "Building with OpenACC (NVHPC)")
    target_compile_definitions(model PRIVATE USE_OPENACC)
    target_compile_options(model PRIVATE
        -fast
        -Minfo=accel
        -Mpreprocess
    )

    # Link flags: need -acc at link time, optionally gpu flags
    target_link_options(model PRIVATE
        -acc
        -gpu=cc80
    )
else()
    message(STATUS "Building with OpenMP-style flags")
    target_compile_definitions(model PRIVATE USE_OPENMP)
    target_compile_options(model PRIVATE
        -Ofast
        -fPIC
        -fopenmp
        -cpp
    )

    # Optional: also use CMake's OpenMP discovery if available
    find_package(OpenMP)
    if(OpenMP_Fortran_FOUND)
        target_link_libraries(model PRIVATE OpenMP::OpenMP_Fortran)
    endif()
endif()

# ------------------------------------------------------------
# Link options / libraries
# ------------------------------------------------------------
# NetCDF Fortran & C libs from nf-config / nc-config
# NF_FLIBS typically contains things like: -L<dir> -lnetcdff -lnetcdf ...
link_directories(${NC_LIBDIR})
target_link_libraries(model PRIVATE ${NF_FLIBS} MPI::MPI_Fortran)

# ------------------------------------------------------------
# Generate documentation
# ------------------------------------------------------------
find_package(Doxygen)
if(Doxygen_FOUND)
    set(DOXYGEN_IN  ${CMAKE_SOURCE_DIR}/doc/Doxyfile) # original config file
    set(DOXYGEN_OUT ${CMAKE_BINARY_DIR}/doc/Doxyfile) # where to store inside 'build' folder
    # set(DOXYGEN_MD_IN  ${CMAKE_SOURCE_DIR}/doc/README.md) # original config file
    # set(DOXYGEN_MD_OUT ${CMAKE_BINARY_DIR}/doc/README.md) # where to store inside 'build' folder
    set(DOXYGEN_STYLE1_IN  ${CMAKE_SOURCE_DIR}/doc/custom.css) # original config file
    set(DOXYGEN_STYLE1_OUT ${CMAKE_BINARY_DIR}/doc/custom.css) # where to store inside 'build' folder
    set(DOXYGEN_STYLE2_IN  ${CMAKE_SOURCE_DIR}/doc/custom_dark_theme.css) # original config file
    set(DOXYGEN_STYLE2_OUT ${CMAKE_BINARY_DIR}/doc/custom_dark_theme.css) # where to store inside 'build' folder

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    # configure_file(${DOXYGEN_MD_IN} ${DOXYGEN_MD_OUT} @ONLY)
    configure_file(${DOXYGEN_STYLE1_IN} ${DOXYGEN_STYLE1_OUT} @ONLY)
    configure_file(${DOXYGEN_STYLE2_IN} ${DOXYGEN_STYLE2_OUT} @ONLY)

    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc # where to save the generated docs
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )
endif()

# ------------------------------------------------------------
# Test setup
# ------------------------------------------------------------
enable_testing()

# If reference NetCDF files are in the source dir, copy them to the build dir
# so the test can find them.
file(GLOB REF_NC_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/output-serial.nc"
    "${CMAKE_CURRENT_SOURCE_DIR}/output-serial-optimized.nc"
)
if(REF_NC_FILES)
    file(COPY ${REF_NC_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

add_test(
    NAME model_test
    COMMAND bash -c
    "echo \"TEST: Running model...\"; \
     OUTPUT=\$(mpirun -np 4 ./model 100 1000 10); \
     MASS=\$(echo \"\${OUTPUT}\" | awk '/Fractional Delta Mass/{print \$NF}'); \
     ENERGY=\$(echo \"\${OUTPUT}\" | awk '/Fractional Delta Energy/{print \$NF}'); \
     awk -v m=\"\${MASS}\" -v e=\"\${ENERGY}\" 'BEGIN { \
       mass = m + 0; energy = e + 0; \
       if (mass < 1e-13 && energy < 1e-3) { \
         print \"TEST PASSED.\"; exit 0; \
       } else { \
         print \"TEST FAILED: thresholds exceeded.\"; \
         printf \"Fractional Delta Mass   = %e  (limit 1e-13)\\n\", mass; \
         printf \"Fractional Delta Energy = %e  (limit 1e-3)\\n\", energy; \
         exit 1; \
       } \
     }'; \
    #  echo \"TEST: Compare output with reference output file...\"; \
    #  IMPROVEMENT as of now we need to manually create virtual env. see `requirements.txt` for more details.
    #  TODO create virtual env automatically, install packages and use it before running the test. will leave this for now.
    #  pyenv/bin/python3 ${CMAKE_CURRENT_SOURCE_DIR}/nccmp3.py output-serial.nc output-serial-optimized.nc output.nc"
)

# Ensure test runs in the build directory (where ./model and output.nc live)
set_tests_properties(model_test PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)