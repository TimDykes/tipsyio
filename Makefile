#-----------------------------------#
#			Demo Makefile			#
#-----------------------------------#

# Compiler
CXX = CC

# Target
TARGET = iotipsy

# Debugging
# Do make DEBUG=1 for debug build, otherwise build with level 2 optimization
ifeq ($(DEBUG), 1)
DBG = -g -O0 
else
DBG = -O2
endif

ifeq ($(MPI),1)
USE_MPI = -DUSE_MPI
endif

# Options and includes
OPTS = -Wall $(DBG) $(USE_MPI)
INCL_PATHS = -I./ 

# Flags for c++ files
CPPFLAGS = $(OPTS) $(INCL_PATHS)  

# Link flags & libraries
LDFLAGS = $(OPTS)
LDLIBS =

# Objects and includes
SRCS = $(wildcard *.cpp) 
OBJS = $(patsubst %.cpp,%.o,$(SRCS))
INCL = $(wildcard *.h) Makefile

.PHONEY: all clean

# Rules
all: $(TARGET) 

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) $(LDLIBS) -o $@

clean: 
	@printf "Cleaning: \n"
	@find . -type f -name '*.o' -print0 | xargs -0 -I % sh -c 'printf "% "; rm -f %'
	rm -f $(TARGET)
