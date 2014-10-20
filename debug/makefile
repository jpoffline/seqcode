# Compilation flags

OPTIMIZE = -O3 -march=native -flto
WARNINGS = -Wall -Wextra -Wno-unused-parameter

CC = g++

CFLAGS = $(OPTIMIZE) $(WARNINGS) -std=gnu++03
INC = -I/usr/local/include
LIBS = -L/usr/local/lib -lm -lgsl -lgslcblas -L/usr/local/boost -lboost_filesystem -lboost_system -lboost_timer

# Make rules

OBJECTS = integrate.o quintessence.o basicdump.o lambdaCDM.o simplecheck.o model.o linearw.o kessence.o process.o KGB.o params.o evolve_sequester.o fr.o inireader.o print2memory.o
OBJECTMAIN = main_new.o 

.PHONY: all
all: main 

main: $(OBJECTS) $(OBJECTMAIN)
	$(CC) $(OBJECTS) $(OBJECTMAIN) $(LIBS) -o ../$@



clean:
	-rm -rf $(OBJECTS) $(OBJECTMAIN) ../main .deps

ifeq ($(USER),jolyon)
# Fancy dependency stuff: stores the dependencies of each source file, so that if headers are updated, it knows which to recompile

DEPDIR = .deps

%.o : ../source/%.cpp
	@mkdir -p .deps && $(CC) -MM $(CFLAGS) $< -MP -MF $(DEPDIR)/$*.P
	$(COMPILE.c) -o $@ -c $(INC) $<

-include $(OBJECTS:%.o=$(DEPDIR)/%.P)
-include $(OBJECTMAIN:%.o=$(DEPDIR)/%.P)
else
%.o: ../source/%.cpp
	$(CC) -c $(INC) $< -o $@
endif
