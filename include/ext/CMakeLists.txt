project(ext-alpscore-headers CXX)

file(GLOB EXT_ALPSCORE_HEADERS alpscore/*.h)

add_custom_target(ext-alpscore-headers SOURCES
    ${EXT_ALPSCORE_HEADERS})
