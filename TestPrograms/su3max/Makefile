
TOPDIR=../..
TESTDIR=$(TOPDIR)/TestPrograms

TESTS = check_element_spinor_transform \
        check_element_gauge_transform \
        check_indexing_transforms \
        check_dirac_operator

LIBS += $(TOPDIR)/LibHR/libhr_su3max.a
CPPFLAGS += -DGAUGE_SUN -DREPR_FUNDAMENTAL -DNG=3

include $(TESTDIR)/MkRules

clean::
	rm -f $(TESTS) *.o
