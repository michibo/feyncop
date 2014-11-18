INSTALL_PATH?=${HOME}

INSTALL_FILES=nauty_wrapper.so geng multig graph.py weighted_graph.py hopf_graph.py phi_k_gen.py phi_34_gen.py qed_gen.py qcd_gen.py combinatorics.py powerseries.py stuff.py nauty_ctrl.py feyncop feyngen


all : nauty_wrapper.so geng multig

nauty_wrapper.so : nauty_wrapper.c ../nauty/nautil.c ../nauty/nauty.c ../nauty/naugraph.c ../nauty/schreier.c ../nauty/naurng.c
	gcc -shared -O2 -lpython2.7 -fPIC $^ -o $@

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
