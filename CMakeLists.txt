cmake_minimum_required(VERSION 3.18)
project(EDLib-alpscore LANGUAGES CXX
        VERSION 2.0)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(CXXStandard)
use_cxx17()

# When EDLib-alpscore is added via FetchContent from EDLib, the EDLib::EDLib
# target already exists in the parent build tree. Otherwise this is a
# standalone build and we need to find the installed EDLib package.
if(NOT TARGET EDLib::EDLib)
    find_package(EDLib REQUIRED)
endif()

find_package(ALPSCore COMPONENTS params hdf5 gf REQUIRED)

add_subdirectory(include/edlib)
add_subdirectory(include/ext)

include(GNUInstallDirs)

add_library(${PROJECT_NAME} INTERFACE)
add_library(EDLib::Alpscore ALIAS ${PROJECT_NAME})

target_include_directories(${PROJECT_NAME}
    INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
              $<BUILD_INTERFACE:${ALPSCore_INCLUDES}>
              $<INSTALL_INTERFACE:include>)

target_link_libraries(${PROJECT_NAME} INTERFACE
    EDLib::EDLib
    ${ALPSCore_LIBRARIES})

target_compile_definitions(${PROJECT_NAME} INTERFACE EDLIB_WITH_ALPSCORE)

# Only run install when this project is top-level. When FetchContent'd by
# EDLib, install rules are driven by the parent project.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    install(TARGETS ${PROJECT_NAME}
            EXPORT ${PROJECT_NAME}_targets
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

    include(CMakePackageConfigHelpers)
    write_basic_package_version_file("${PROJECT_NAME}ConfigVersion.cmake"
                                     VERSION ${PROJECT_VERSION}
                                     COMPATIBILITY SameMajorVersion)
    configure_package_config_file(
        "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
        "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
        INSTALL_DESTINATION
        ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)
    install(EXPORT ${PROJECT_NAME}_targets
            FILE ${PROJECT_NAME}Targets.cmake
            NAMESPACE EDLib::
            DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)
    install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
                  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
            DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

    install(DIRECTORY include/edlib/alpscore DESTINATION include/edlib
            FILES_MATCHING PATTERN "*.h")
    install(DIRECTORY include/ext/alpscore   DESTINATION include/ext
            FILES_MATCHING PATTERN "*.h")
endif()
