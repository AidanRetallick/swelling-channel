# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
noinst_PROGRAMS = rectangular_plate_koiter_steigmann           \
                  testint_rectangular_plate_koiter_steigmann   \
                  inflate_rectangular_plate_koiter_steigmann   \
		  perturb_rectangular_plate_koiter_steigmann   \
                  prestrain_rectangular_plate_koiter_steigmann \
		  symmetric_rectangular_plate_koiter_steigmann \
                  rectangular_pressure_cycle_koiter_steigmann  \
                  infinite_strip_koiter_steigmann              \
                  euler_column_buckling_koiter_steigmann
include $(top_srcdir)/config/makefile_templates/demo_drivers
#---------------------------------------------------------------------------

# Local sources that each code depends on:
rectangular_plate_koiter_steigmann_SOURCES = \
 rectangular_plate_koiter_steigmann.cc

inflate_rectangular_plate_koiter_steigmann_SOURCES = \
 inflate_rectangular_plate_koiter_steigmann.cc

perturb_rectangular_plate_koiter_steigmann_SOURCES = \
 perturb_rectangular_plate_koiter_steigmann.cc

prestrain_rectangular_plate_koiter_steigmann_SOURCES = \
 prestrain_rectangular_plate_koiter_steigmann.cc

symmetric_rectangular_plate_koiter_steigmann_SOURCES = \
 symmetric_rectangular_plate_koiter_steigmann.cc

rectangular_pressure_cycle_koiter_steigmann_SOURCES = \
 rectangular_pressure_cycle_koiter_steigmann.cc

infinite_strip_koiter_steigmann_SOURCES = \
 infinite_strip_koiter_steigmann.cc

euler_column_buckling_koiter_steigmann_SOURCES = \
 euler_column_buckling_koiter_steigmann.cc
#---------------------------------------------------------------------------

rectangular_plate_koiter_steigmann_LDADD = -L@libdir@ -lc1_koiter_steigmann \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
inflate_rectangular_plate_koiter_steigmann_LDADD = -L@libdir@ -lc1_koiter_steigmann \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
perturb_rectangular_plate_koiter_steigmann_LDADD = -L@libdir@ -lc1_koiter_steigmann \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
prestrain_rectangular_plate_koiter_steigmann_LDADD = -L@libdir@ -lc1_koiter_steigmann \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
symmetric_rectangular_plate_koiter_steigmann_LDADD = -L@libdir@ -lc1_koiter_steigmann \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
rectangular_pressure_cycle_koiter_steigmann_LDADD = -L@libdir@ -lc1_koiter_steigmann \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
infinite_strip_koiter_steigmann_LDADD = -L@libdir@ -lc1_koiter_steigmann \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
euler_column_buckling_koiter_steigmann_LDADD = -L@libdir@ -lc1_koiter_steigmann \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Local sources that Jack's own code depends on: This code also uses
# objects from Jack's own library. The source code for this library
# lives in user_src/jack_cool and gets compiled in to the jack_cool 
# library. If we use it, we have to link against it...
#jacks_own_code_SOURCES = jacks_own_code.cc

# Required libraries: Jack's own library ("jack_cool") and oomph-lib's 
# and "generic" and "poisson" libraries. Following the logic
# explained above, Jack's own library should go first since it may
# (or may not) depend on objects in the other libraries, whereas they will
# certainly not depend on objects in his. 
# jacks_own_code_LDADD = -L@libdir@ -ljack_cool -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------


# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.

AM_CPPFLAGS += -I@includedir@ -Wall -ggdb3 -O3
#EXTRA_DIST = Circle1.1.ele Circle1.1.node Circle1.1.poly
