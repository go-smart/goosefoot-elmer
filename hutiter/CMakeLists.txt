
cONFIGURE_FILE(config.h.cmake config.h)

INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}")
INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}/hutiter")
INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}/hutiter/src")

ADD_SUBDIRECTORY(src)
