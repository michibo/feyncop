
# Standard things

sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

PACKAGE_BUILT_$(d) := $(d)/../feyncop_built.tar.gz
PACKAGE_$(d) := $(d)/../feyncop.tar.gz

TGTS_$(d)   := $(PACKAGE_BUILT_$(d)) $(PACKAGE_$(d))

pack : $(TGTS_$(d))

SRCS_$(d)   := $(addprefix $(d)/, graph.py weighted_graph.py hopf_graph.py phi_k_gen.py phi_34_gen.py qed_gen.py qcd_gen.py combinatorics.py powerseries.py stuff.py nauty_ctrl.py feyncop feyngen)

EXTRA_SRCS_$(d) := $(addprefix $(d)/, nauty_wrapper.c Makefile)

MISC_$(d)   := $(addprefix $(d)/, README)

BINARIES_$(d)   := $(addprefix $(d)/, nauty_wrapper.so geng multig)

$(PACKAGE_BUILT_$(d)) : $(SRCS_$(d)) $(MISC_$(d)) $(BINARIES_$(d))
		tar czf $@ $^

$(PACKAGE_$(d)) : $(SRCS_$(d)) $(MISC_$(d)) $(EXTRA_SRCS_$(d))
		tar czf $@ $^

CLEAN		:= $(CLEAN) $(TGTS_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
