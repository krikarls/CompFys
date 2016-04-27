# the compiler
CC = gcc
CXX = g++


# forces removal without prompt
RM = rm -f

# compiler flags:
# -g debugging information to executable
# -Wall most compiler warnings enabled
# -Werror turns warnings into errors
# -O3 optimiation level 3
# -pedantic warns about GCC spesific options
# -MMD & -MP makes automatic dependancy files
CPPFLAGS = -g -std=c++11 -Wall -Wextra -Werror -O3 -pedantic -MMD -MP -fopenmp

# linker flags
LDFLAGS = -g -Wall


# libraries to link into executable
LDLIBS =  -fopenmp

# build executable
TARGET = project4
CPP_FILES = main.cpp lattice.cpp ising.cpp

SRCS = $(addprefix src/,$(notdir $(CPP_FILES)))
OBJS = $(patsubst src/%.cpp,obj/%.o,$(SRCS))


TEST_TARGET = testit
TEST_FILES = tests.cpp catch.cpp lattice.cpp ising.cpp
TEST_SRCS = $(addprefix src/,$(notdir $(TEST_FILES)))
TEST_OBJS = $(patsubst src/%.cpp,obj/%.o,$(TEST_SRCS))

ANIMATE_TARGET = animation
ANIMATE_FILES = animate.cpp lattice.cpp ising.cpp 
ANIMATE_SRCS = $(addprefix src/,$(notdir $(ANIMATE_FILES)))
ANIMATE_OBJS = $(patsubst src/%.cpp,obj/%.o,$(ANIMATE_SRCS))
ANIMATE_LIBS = -lX11 -lpthread


default: $(TARGET) $(TEST_TARGET)
	./testit # lazy auto-testing

all: $(TARGET) $(TEST_TARGET) $(ANIMATE_TARGET) docs


$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

obj/%.o: src/%.cpp
	@mkdir -p obj
	$(CXX) $(CPPFLAGS) -c -o $@ $<
	


test: $(TEST_TARGET)


$(TEST_TARGET): $(TEST_OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)


animate: $(ANIMATE_TARGET)

$(ANIMATE_TARGET): $(ANIMATE_OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(ANIMATE_LIBS)
	

docs:
	doxygen Doxyfile

.phony: all clean dist-clean test animation docs

clean:
	$(RM) $(patsubst %.o,%.d,$(OBJS))
	$(RM) $(OBJS)
	$(RM) $(patsubst %.o,%.d,$(TEST_OBJS))
	$(RM) $(TEST_OBJS)
	$(RM) $(patsubst %.o,%.d,$(ANIMATE_OBJS))
	$(RM) $(ANIMATE_OBJS)

dist-clean: clean
	$(RM) $(TARGET)
	$(RM) $(TEST_TARGET)
	$(RM) $(ANIMATE_TARGET)
	@rmdir --ignore-fail-on-non-empty obj
	$(RM) -r html




