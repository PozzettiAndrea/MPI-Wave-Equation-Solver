.PHONY: main clean

CXX = mpicxx
CXXFLAGS = -std=c++17
OPTFLAGS = -O3 -march=native

main: main.cpp Funcs.cpp
	$(CXX)  $(OPTFLAGS) -o $@ main.cpp Funcs.cpp Funcs.h $(CXXFLAGS)

clean:
	-rm -rf main