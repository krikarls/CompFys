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
LDFLAGS = -g -Wall


# libraries to link into executable
LDLIBS = -larmadillo

# build executable
TARGET = project3
SRCS = main.cpp gaussian_quad.cpp monte_carlo.cpp
OBJS = $(subst .cpp,.o,$(SRCS))


all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(TARGET) $(OBJS) $(LDLIBS)



clean:
	$(RM) $(OBJS)

