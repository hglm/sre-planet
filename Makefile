
CONFIGURATION = DEBUG

# CFLAGS with optional tuning for CPU
OPTCFLAGS = -Ofast -ffast-math
ifeq ($(TARGET_CPU), CORTEX_A7)
OPTCFLAGS += -mcpu=cortex-a7 -mfpu=vfpv4
endif
ifeq ($(TARGET_CPU), CORTEX_A8)
OPTCFLAGS += -mcpu=cortex-a8
endif
ifeq ($(TARGET_CPU), CORTEX_A9)
OPTCFLAGS += -mcpu=cortex-a9
endif
ifeq ($(CONFIGURATION), DEBUG)
OPTCFLAGS = -ggdb
endif
CFLAGS = -Wall -pipe $(OPTCFLAGS)

# pkgconfig dependencies of SRE will be automatically generated.
PKG_CONFIG_CFLAGS = `pkg-config --cflags sre`
PKG_CONFIG_LIBS = `pkg-config --libs sre`

MODULE_OBJECTS = sre-planet.o mesh.o
PROGRAM = sre-planet

$(PROGRAM) : $(MODULE_OBJECTS)
	g++ $(MODULE_OBJECTS) -o $(PROGRAM) $(PKG_CONFIG_LIBS)

.cpp.o :
	$(CC) -c $(CFLAGS) $(PKG_CONFIG_CFLAGS) $< -o $@

clean :
	rm -f $(MODULE_OBJECTS)
	rm -f sre-planet

dep :
	rm -f .depend
	make .depend

.depend: Makefile
	echo '# Module dependencies' >> .depend
	g++ -MM $(patsubst %.o,%.cpp,$(MODULE_OBJECTS)) >> .depend

include .depend
