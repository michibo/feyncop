
# Standard things

sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

TGTS_$(d)   := $(d)/../feyncop_built.tar.gz

pack : $(TGTS_$(d))

SRCS_$(d)   := $(addprefix $(d)/, graph.py weighted_graph.py hopf_graph.py phi_k_gen.py phi_34_gen.py qed_gen.py qcd_gen.py combinatorics.py powerseries.py stuff.py nauty_ctrl.py feyncop feyngen)

MISC_$(d)   := $(addprefix $(d)/, nauty_wrapper.so geng multig README)

$(TGTS_$(d)) : $(SRCS_$(d)) $(MISC_$(d))
		tar czf $@ $^

DEPS_$(d)	:= $(TGTS_$(d):%.pdf=%.d)
SRCS_$(d)   := $(TGTS_$(d):%.pdf=%.tex)

CLEAN		:= $(CLEAN) $(TGTS_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
