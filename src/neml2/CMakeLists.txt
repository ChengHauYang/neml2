# ----------------------------------------------------------------------------
# Helper function to add a submodule
# ----------------------------------------------------------------------------
function(neml2_add_submodule mname TYPE mdir)
  file(GLOB_RECURSE msrcs CONFIGURE_DEPENDS ${mdir}/*.cxx)
  file(GLOB_RECURSE mheaders CONFIGURE_DEPENDS
          ${mdir}/include/**/*.h
          ${NEML2_SOURCE_DIR}/include/neml2/${mdir}/*.h
  )
  add_library(${mname} ${TYPE} ${msrcs})

  target_sources(${mname}
    PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${NEML2_SOURCE_DIR}/include
    FILES
    ${mheaders}
  )

  set_target_properties(${mname} PROPERTIES OUTPUT_NAME "neml2_${mname}$<IF:$<CONFIG:Release>,,_$<CONFIG>>")

  if(NEML2_WHEEL)
    set_target_properties(${mname} PROPERTIES INSTALL_RPATH "${INSTALL_REL_PATH};${INSTALL_REL_PATH}/../../torch/lib")
  else()
    set_target_properties(${mname} PROPERTIES INSTALL_RPATH "${INSTALL_REL_PATH};${torch_LINK_DIR}")
  endif()

  target_compile_options(${mname} PRIVATE -Wall -Wextra -pedantic)

  # Add library to the interface neml2 target
  target_link_libraries(neml2 INTERFACE ${mname})

  install(TARGETS ${mname}
    EXPORT neml2targets
    LIBRARY COMPONENT libneml2
    FILE_SET HEADERS COMPONENT libneml2
  )

  # Header search paths for all submodules
  target_include_directories(${mname}
    PUBLIC
      $<BUILD_INTERFACE:${NEML2_SOURCE_DIR}/include/neml2>
      $<BUILD_INTERFACE:${NEML2_SOURCE_DIR}/tensors/include>
      $<INSTALL_INTERFACE:include>
  )

endfunction()

# ----------------------------------------------------------------------------
# Submodules
# ----------------------------------------------------------------------------
# libneml2
#
# This library is a wrapper of all neml2 library components.
#
# Since neml2 relies on dynamic loading of object symbols, targets who link
# against neml2 should preferrably use LINK_WHAT_YOU_USE. Otherwise, dlopen-like
# utilities shall be used to load the runtime objects.
add_library(neml2 INTERFACE)

# ----------------------------------------------------------------------------
# Set include paths
# ----------------------------------------------------------------------------
# Use generator expressions ($<...>) here to fix export errors
target_include_directories(neml2 INTERFACE
  # 1. Let everyone find the tensors and misc headers (currently under tensors/include)
  #    Assuming your layout is: path/to/tensors/include/misc/errors.h
  $<BUILD_INTERFACE:${NEML2_SOURCE_DIR}/tensors/include>

  # 2. Make includes like base/OptionSet.h work
  #    Since Registry.h does #include "base/OptionSet.h"
  #    we need to add .../include/neml2 to the search paths
  $<BUILD_INTERFACE:${NEML2_SOURCE_DIR}/include/neml2>
)

if(NEML2_JSON)
  target_compile_definitions(misc PUBLIC NEML2_HAS_JSON)
  target_link_libraries(misc PUBLIC nlohmann_json::nlohmann_json)
endif()

if(NEML2_CSV)
  target_compile_definitions(misc PUBLIC NEML2_HAS_CSV)
  target_link_libraries(misc PUBLIC csv)
endif()

if(NEML2_PCH)
target_precompile_headers(tensors PUBLIC
        <misc/types.h>
        <misc/assertions.h>
)
endif()

# libneml2_base
neml2_add_submodule(base SHARED base)
target_link_libraries(base PUBLIC tensors hit)

if(NEML2_PCH)
  target_precompile_headers(base PUBLIC
    <neml2/base/Registry.h>
    <neml2/base/Factory.h>
    <neml2/base/Option.h>
    <neml2/base/TensorName.h>
  )
endif()

# libneml2_user_tensor
neml2_add_submodule(user_tensor SHARED user_tensors)
target_link_libraries(user_tensor PUBLIC tensors base)
target_link_options(user_tensor INTERFACE ${CMAKE_CXX_LINK_WHAT_YOU_USE_FLAG})

# libneml2_solver
neml2_add_submodule(solver SHARED solvers)
target_link_libraries(solver PUBLIC tensors base)
target_link_options(solver INTERFACE ${CMAKE_CXX_LINK_WHAT_YOU_USE_FLAG})

# libneml2_model
neml2_add_submodule(model SHARED models)
target_link_libraries(model PUBLIC solver tensors base)
target_link_options(model INTERFACE ${CMAKE_CXX_LINK_WHAT_YOU_USE_FLAG})

if(NEML2_PCH)
  target_precompile_headers(model PUBLIC
    <neml2/models/Model.h>
  )
endif()

# libneml2_driver
neml2_add_submodule(driver SHARED drivers)
target_link_libraries(driver PUBLIC model tensors base)
target_link_options(driver INTERFACE ${CMAKE_CXX_LINK_WHAT_YOU_USE_FLAG})

# libneml2_dispatcher
if(NEML2_WORK_DISPATCHER)
  neml2_add_submodule(dispatcher SHARED dispatchers)
  target_link_libraries(dispatcher PUBLIC tensors model)
  target_link_libraries(dispatcher PUBLIC timpi)
  target_compile_definitions(dispatcher PUBLIC NEML2_HAS_DISPATCHER)
  target_link_libraries(driver PUBLIC dispatcher)
endif()

# ----------------------------------------------------------------------------
# Version number and hash
# ----------------------------------------------------------------------------
find_package(Git)

if(Git_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --abbrev=0
    WORKING_DIRECTORY ${NEML2_SOURCE_DIR}
    OUTPUT_VARIABLE NEML2_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  file(WRITE ${NEML2_BINARY_DIR}/version ${NEML2_VERSION})

  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
    WORKING_DIRECTORY ${NEML2_SOURCE_DIR}
    OUTPUT_VARIABLE NEML2_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  file(WRITE ${NEML2_BINARY_DIR}/hash ${NEML2_HASH})

  install(FILES
    ${NEML2_BINARY_DIR}/version
    ${NEML2_BINARY_DIR}/hash
    DESTINATION .
    COMPONENT libneml2
  )
endif()

# ----------------------------------------------------------------------------
# Export targets
# ----------------------------------------------------------------------------
install(TARGETS neml2 EXPORT neml2targets LIBRARY COMPONENT libneml2)
install(EXPORT neml2targets NAMESPACE neml2:: DESTINATION share/cmake/neml2 COMPONENT libneml2)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${NEML2_SOURCE_DIR}/cmake/neml2Config.cmake.in
  ${NEML2_BINARY_DIR}/neml2Config.cmake
  INSTALL_DESTINATION share/cmake/neml2
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)
write_basic_package_version_file(
  ${NEML2_BINARY_DIR}/neml2ConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(
  FILES
  ${NEML2_BINARY_DIR}/neml2Config.cmake
  ${NEML2_BINARY_DIR}/neml2ConfigVersion.cmake
  DESTINATION share/cmake/neml2
  COMPONENT libneml2
)
