main: gx_core.h $(addsuffix .so,$(basename $(wildcard *.c)))
#	echo $(addsuffix .so,$(basename $(wildcard *.c)))
.PHONY : main
%.so:%.c gx_core.h
	-rm -f $(addsuffix .o,$(basename $<))
	-rm -f $(addsuffix .so,$(basename $<))
	R CMD SHLIB $<
	-rm -f $(addsuffix .o,$(basename $<))
	

	