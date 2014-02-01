
all : nauty_wrapper.so geng multig

nauty_wrapper.so : nauty_wrapper.c ../nauty/nautil.c ../nauty/nauty.c ../nauty/naugraph.c ../nauty/schreier.c ../nauty/naurng.c
	gcc -shared -O2 -lpython2.7 -fPIC $^ -o $@

geng : ../nauty/geng
	cp $^ $@

multig : ../nauty/multig
	cp $^ $@

clean:
	rm -f multig geng nauty_wrapper.so *.pyc
