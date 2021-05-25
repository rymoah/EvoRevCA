CXX = g++
WARN = -Wall
CFLAGS = $(WARN)
LDFLAGS = $(WARN) -lecf

OBJS = EvalOp.o \
       Primitives.o \
       main.o

BINFILE = revca

.PHONY = all clean

all: $(BINFILE)

$(BINFILE): $(OBJS)
	$(CXX) $(OBJS) $(LDFLAGS) -o $(BINFILE)

%.o: %.cpp %.h
	$(CXX) $(CFLAGS) $< -c -o $@

clean:
	rm -f *~ $(OBJS) $(BINFILE)
