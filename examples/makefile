.PHONY: all rel clean
TARGETS := test_prog example
CFLAGS:= -g -O2 -Wall -Wextra -Wno-unused-function -I../src
all: rel
dbg: CFLAGS := -DNSP_DEBUG -g3 -O0 -fsanitize=address,undefined -Wall -Wextra -Wno-unused-function -I../src
dbg: $(TARGETS)
rel: $(TARGETS)
clean:
	rm $(TARGETS)
