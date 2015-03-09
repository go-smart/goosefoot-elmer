
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(writetest_SRCS writetest.f90)

if(KDE4_BUILD_TESTS)

kde4_add_executable(writetest ${writetest_SRCS})

target_link_libraries(writetest ${KDE4_KDECORE_LIBS} binio)

endif(KDE4_BUILD_TESTS)

########### next target ###############

set(readtest_SRCS readtest.f90)

if(KDE4_BUILD_TESTS)

kde4_add_executable(readtest ${readtest_SRCS})

target_link_libraries(readtest ${KDE4_KDECORE_LIBS} binio)

endif(KDE4_BUILD_TESTS)

########### install files ###############




#original Makefile.am contents follow:

#TESTS = run-writetest run-readtest
#check_PROGRAMS = writetest readtest
#
#run-writetest: writetest.sh
#	cp writetest.sh run-writetest
#	chmod a+x run-writetest
#
#run-readtest: readtest.sh
#	cp readtest.sh run-readtest
#	chmod a+x run-readtest
#
#AM_FCFLAGS = @INCLUDE_MODULE_FLAG@..
#
#writetest_SOURCES = writetest.f90
#writetest_LDADD = ../libbinio.a
#
#readtest_SOURCES = readtest.f90
#readtest_LDADD = ../libbinio.a
#
#EXTRA_DIST = writetest.sh readtest.sh writetest.checkL writetest.checkB
#
#clean-local:
#	$(RM)  writetest.out
