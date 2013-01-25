main : $(addsuffix .so,$(basename $(wildcard *.c)))
.PHONY : main
%.so : %.c
	-rm -f $(addsuffix .o,$(basename $<))
	-rm -f $(addsuffix .so,$(basename $<))
	-R CMD SHLIB $<
	-rm -f $(addsuffix .o,$(basename $<))
