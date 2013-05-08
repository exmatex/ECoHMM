ERL=erl
ERLC=erlc
ERLCFLAGS+=-W +debug_info
ERLS=comd_sup.erl comd_srv.erl hmm.erl hmm_sup.erl hmm_srv.erl
BEAMS=$(ERLS:.erl=.beam)

.PHONY: clean
.SUFFIXES: .beam .erl 

all: $(BEAMS)

.erl.beam:
	$(ERLC) $(ERLCFLAGS) $<

clean:
	rm -f $(BEAMS) 

run:
	$(ERL) 