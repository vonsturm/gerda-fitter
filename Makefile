# Makefile
#
# Author: Luigi Pertoldi - pertoldi@pd.infn.it
# Created: Mon 6 May 2019

FLAGS = -Wall -Wextra -std=c++11 -O3 \
        $$(bat-config --cflags)
LIBS  = $$(bat-config --libs)

gerda-fitter : src/gerda-fitter.cc src/GerdaFitter.cc
	$(CXX) $(FLAGS) -o $@ $^ $(LIBS)

install : gerda-fitter
	install -d $(PREFIX)/bin
	install $< $(PREFIX)/bin
