# test level Makefile.am

AM_CPPFLAGS += -I$(top_srcdir)/src/include
AM_LDFLAGS += -L${abs_top_builddir}/src/lib

LDADD =
LDADD += -lH5hut

noinst_PROGRAMS =

if ENABLE_C
noinst_PROGRAMS += h5u_test h5b_test h5_attach_test

TESTS_ENVIRONMENT = env LD_LIBRARY_PATH=@HDF5_PREFIX@/lib:$(LD_LIBRARY_PATH)

EXTRA_PROGRAMS = h5u_test h5b_test h5_attach_test

h5u_test_SOURCES =	\
		h5u_test.c \
		h5u_write.c \
		h5u_read.c \
		testframe.c \
		testframe.h \
		params.h

h5b_test_SOURCES =	\
		h5b_test.c \
		h5b_write.c \
		h5b_read.c \
		testframe.c \
		testframe.h \
		params.h

h5_attach_test_SOURCES = \
		h5_attach_test.c
endif

clean-local:
	${RM} -f *.h5

