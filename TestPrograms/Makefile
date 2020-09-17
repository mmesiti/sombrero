TOPDIR=..

SUBDIRS := $(shell find * -maxdepth 0 -type d)

LIBS += $(TOPDIR)/LibHR/libhr_su3max.a

tests: $(LIBS)
	touch $(TOPDIR)/LibHR/.lib.updated
	for i in $(SUBDIRS) ; \
	do \
	  $(MAKE) -C $$i ;\
	done
	rm -f touch $(TOPDIR)/LibHR/.lib.updated
