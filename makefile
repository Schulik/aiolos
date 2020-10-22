CXX = g++ -std=c++17
CPPFLAGS =      # put pre-processor settings (-I, -D, etc) here
CXXFLAGS = -Wall -Wextra #-pg or -g put compiler settings here
LDFLAGS = -lm -flto      # put linker settings here
BFLAGS = -I. -O3 -g -flto #-pg or -g

PROBLEM=default

##SRC = main.cpp advection.cpp source.cpp
SRC = $(wildcard *.cpp)
SRC += problems/$(PROBLEM).cpp
OBJ = $(SRC:.cpp=.o)
INC = $(SRC:.cpp=.h)


TEST_OBJ = $(subst main.o, test_files/main.o, $(OBJ))



#euler_hllc: $(OBJ)
#	$(CXX) -o $@ $(OBJ) $(CXXFLAGS) $(LDFLAGS) $(BFLAGS)

aiolos: $(OBJ) makefile
	$(CXX) -o $@ $(OBJ) $(CXXFLAGS) $(LDFLAGS)
	#./$@

%.o: %.cpp makefile aiolos.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< $(BFLAGS) -o $@ 

.PHONY: tests

tests: $(TEST_OBJ) makefile aiolos.h
	$(CXX) -o  $@  $(TEST_OBJ) $(CXXFLAGS) $(LDFLAGS)
	./tests > /dev/null
	cd test_files ; python test_shock_tube.py
	cd test_files ; python test_steady_state.py

clean:
	rm -f *.o test_files/*.o problems/*.o
