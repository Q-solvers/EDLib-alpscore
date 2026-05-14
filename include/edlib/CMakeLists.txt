project(edlib-alpscore-headers CXX)

file(GLOB EDLIB_ALPSCORE_HEADERS alpscore/*.h)

add_custom_target(edlib-alpscore-headers SOURCES
    ${EDLIB_ALPSCORE_HEADERS})
