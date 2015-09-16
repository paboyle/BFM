SUBDIRS=

if BUILD_QDP
SUBDIRS+= \
	wilson \
        clover \
        clover_threaded \
        clover_cg \
	wilson_vec \
	linalg \
	dwf_rb5d_dperp \
	dwf_rb5d \
	dwf_prec \
	dwf_mooee \
	dwf_mooeeinv \
	dwf_rb5d_cg \
	dwf_cg \
	dwf_rb5d_prop \
	wilson_cg
else
SUBDIRS+=
	wilson_nochroma \
	wilson_singlenode \
        dwf_singlenode \
	dwf_nochroma \
	dwf_cg_nochroma 
endif


.PHONY: $(SUBDIRS)

recurse: $(SUBDIRS)

install:

all: recurse 

clean: recurse 

distclean: recurse 


$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)


