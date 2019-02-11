# Makefile du programme C++ Analyse

SHELL = /bin/sh
prefix = /usr/local
bindir = $(prefix)/bin
mandir = $(prefix)/man/man1
EXEC = Analyse
SH = lancer_analyse.sh
CXXFLAGS = -O2 -std=c++0x -march=native
LIBS = -lz -pthread

all: $(EXEC) $(SH) $(SH).1

$(EXEC): $(EXEC).cpp eval.o eval.h
	g++ $(EXEC).cpp eval.o $(CXXFLAGS) $(LIBS) -o $@

eval.o: eval.cpp eval.h
	g++ $< $(CXXFLAGS) -c -o $@

install: all
	install $(EXEC) $(bindir)/$(EXEC)
	install $(SH) $(bindir)/$(SH)
	@mkdir -p $(mandir)
	install $(SH).1 $(mandir)/$(SH).1
	@ln -f -s $(mandir)/$(SH).1 $(mandir)/$(EXEC).1

clean:
	@rm -f *.o $(EXEC)

remove: 
	@rm -f $(bindir)/$(SH) $(bindir)/$(EXEC) \
	$(mandir)/$(SH).1 $(mandir)/$(EXEC).1