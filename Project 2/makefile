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
CPPFLAGS = -g -std=c++11 -Wall -Wextra -Werror -O3 -pedantic -MMD -MP

# linker flags
LDFLAGS = -g -Wall


# libraries to link into executable
LDLIBS = -larmadillo -llapack -lopenblas -lunittest++

# build executable
TARGET = project2
CPP_FILES = project2.cpp Jacobi_rotation.cpp test.cpp filewriter.cpp helper_files.cpp


SRCS = $(addprefix src/,$(notdir $(CPP_FILES)))


OBJS = $(patsubst src/%.cpp,obj/%.o,$(SRCS))



all: $(TARGET) 


$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

obj/%.o: src/%.cpp
	@mkdir -p obj
	$(CXX) $(CPPFLAGS) -c -o $@ $<
	

.phony: all clean dist-clean

clean:
	$(RM) $(patsubst %.o,%.d,$(OBJS))
	$(RM) $(OBJS)

dist-clean: clean
	$(RM) $(TARGET)
	@rmdir --ignore-fail-on-non-empty obj 




