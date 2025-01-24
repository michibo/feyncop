INSTALL_PATH?=${HOME}

INSTALL_FILES=nauty_wrapper.so geng multig graph.py weighted_graph.py hopf_graph.py phi_k_gen.py phi_34_gen.py qed_gen.py qcd_gen.py combinatorics.py powerseries.py stuff.py nauty_ctrl.pyx parsefg.py outputfg.py feyncop feyngen

PY_CFLAGS  := $(shell python-config --cflags)


all : nauty_wrapper.so nauty_ctrl.c geng multig

nauty_wrapper.so : nauty_wrapper.c ../nauty/nautil.c ../nauty/nauty.c ../nauty/naugraph.c ../nauty/schreier.c ../nauty/naurng.c
	gcc -shared -O2 ${PY_CFLAGS} -fPIC $^ -o $@

nauty_ctrl.c : nauty_ctrl.pyx
	cythonize -a -i -3 $<

geng : ../nauty/geng
	cp $^ $@

multig : ../nauty/multig
	cp $^ $@

install : ${INSTALL_FILES}
	rm -f $(addsuffix c, $(filter %.py, ${INSTALL_FILES}))
	cp ${INSTALL_FILES} ${INSTALL_PATH}/bin

.PHONY: uninstall clean
uninstall:
	rm -f $(addsuffix c, $(filter %.py, ${INSTALL_FILES}))
	rm -f $(addprefix ${INSTALL_PATH}/bin/, ${INSTALL_FILES})

clean:
	rm -f multig geng nauty_wrapper.so *.pyc
