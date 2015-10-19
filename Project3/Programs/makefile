# the compiler
CC = gcc
CXX = g++

# forces removal without prompt
RM = rm -f

# compiler flags:
# -g debugging information to executable
# -Wall most compiler warnings enabled
CPPFLAGS = -g -std=c++11 -Wall -Wextra -Werror -O3 -pedantic -fopenmp
# linker flags
LDFLAGS = -g -Wall -fopenmp


# libraries to link into executable
LDLIBS = -larmadillo

# build executable
TARGET = project3 
SRCS = main.cpp gaussian_quadrature.cpp monte_carlo_parallell.cpp
OBJS = $(subst .cpp,.o,$(SRCS))


all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o  $@ $^ $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) -c -o $@ $<



clean:
	$(RM) $(OBJS)

