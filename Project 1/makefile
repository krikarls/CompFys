# the compiler
CC = gcc
CXX = g++

# forces removal without prompt
RM = rm -f

# compiler flags:
# -g debugging information to executable
# -Wall most compiler warnings enabled
CPPFLAGS = -g -std=c++11 -Wall -Wextra -Werror -O3 -pedantic
# linker flags
LDFLAGS = -g


# libraries to link into executable
LDLIBS = -larmadillo -lsuperlu -llapack -lopenblas

# build executable
TARGET = project1
SRCS = $(TARGET).cpp tridiag_solver.cpp
OBJS = $(subst .cpp,.o,$(SRCS))


all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(TARGET) $(OBJS) $(LDLIBS)

depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS)

dist-clean: clean
	$(RM) *~.depend

include .depend
