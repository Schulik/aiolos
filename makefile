

ifeq ($(SYSTEM), intel)
CXX = icpc -std=c++17
CXXFLAGS = -Wall -Wextra -xHost -ipo #-pg or -g put compiler settings here
BFLAGS = -I. -O3 -g #or -g
else 
ifeq ($(SYSTEM), clang)
CXX = clang++ -std=c++17 
CXXFLAGS = -Wall -Wextra -march=native -flto  #-pg or -g put compiler settings here
BFLAGS = -I. -O3 -g #or -g
else
CXX = g++ -std=c++17 
CXXFLAGS = -Wall -Wextra -march=native -flto  #-pg or -g put compiler settings here
BFLAGS = -I. -O3 #or -g
endif
endif
CPPFLAGS = -I/usr/include/eigen3 -DNDEBUG # put pre-processor settings (-I, -D, etc) here
LDFLAGS = -lgsl -lgslcblas -lm # put linker settings here


PROBLEM=default
NUM_SPECIES=Eigen::Dynamic

CPPFLAGS += -DNUM_SPECIES=$(NUM_SPECIES)

##SRC = main.cpp advection.cpp source.cpp
SRC = $(wildcard *.cpp)
SRC += problems/$(PROBLEM).cpp
OBJ = $(SRC:.cpp=.o)
INC = $(SRC:.cpp=.h)


TEST_OBJ = $(subst main.o, test_files/main.o, $(OBJ))



aiolos: $(OBJ) makefile
	$(CXX) -o $@ $(OBJ) $(CXXFLAGS) $(LDFLAGS)

%.o: %.cpp makefile aiolos.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< $(BFLAGS) -o $@ 

.PHONY: tests

tests: $(TEST_OBJ) makefile aiolos.h
	$(CXX) -o  $@  $(TEST_OBJ) $(CXXFLAGS) $(LDFLAGS)
	./tests > /dev/null
	cd test_files ; python3 test_shock_tube.py -p
	cd test_files ; python3 test_steady_state.py
	cd test_files ; python3 test_drag.py
	cd test_files ; python3 test_dustywave.py -p 
	cd test_files ; python3 test_dustyshock.py -p
	cd test_files ; python3 test_conservation.py
	cd test_files ; python3 test_irradiation.py -p
	cd test_files ; python3 test_coll_heating.py

clean:
	rm -f *.o test_files/*.o problems/*.o test_files/*dat
