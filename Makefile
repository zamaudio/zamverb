#!/usr/bin/make -f

PREFIX ?= /usr
LIBDIR ?= lib
LV2DIR ?= $(PREFIX)/$(LIBDIR)/lv2

#OPTIMIZATIONS ?= -msse -msse2 -mfpmath=sse -ffast-math -fomit-frame-pointer -O3 -fno-finite-math-only
OPTIMIZATIONS ?= -ffast-math -fomit-frame-pointer -O3 -fno-finite-math-only

#OPTIMIZATIONS ?=

LDFLAGS ?= -Wl,--as-needed -lm
CXXFLAGS ?= $(OPTIMIZATIONS) -Wall -g
CFLAGS ?= $(OPTIMIZATIONS) -Wall -g

###############################################################################
BUNDLE = zamverb.lv2

CC = gcc
CXX = g++
CXXFLAGS += -fPIC -DPIC
CFLAGS += -fPIC -DPIC

UNAME=$(shell uname)
ifeq ($(UNAME),Darwin)
  LIB_EXT=.dylib
  LDFLAGS += -dynamiclib
else
  LDFLAGS += -shared -Wl,-Bstatic -Wl,-Bdynamic
  LIB_EXT=.so
endif


ifeq ($(shell pkg-config --exists lv2 lv2-plugin || echo no), no)
  $(error "LV2 SDK was not found")
else
  LV2FLAGS=`pkg-config --cflags --libs lv2 lv2-plugin`
endif

$(BUNDLE): manifest.ttl zamverb.ttl zamverb$(LIB_EXT)
#zamverb_gui$(LIB_EXT)
	rm -rf $(BUNDLE)
	mkdir $(BUNDLE)
	cp manifest.ttl zamverb.ttl zamverb$(LIB_EXT) $(BUNDLE)

zamverb$(LIB_EXT): zamverb.c
	$(CC) -o zamverb$(LIB_EXT) \
		$(CXXFLAGS) \
		zamverb.c \
		$(LV2FLAGS) $(LDFLAGS)

install: $(BUNDLE)
	install -d $(DESTDIR)$(LV2DIR)/$(BUNDLE)
	install -t $(DESTDIR)$(LV2DIR)/$(BUNDLE) $(BUNDLE)/*

uninstall:
	rm -rf $(DESTDIR)$(LV2DIR)/$(BUNDLE)

clean:
	rm -rf $(BUNDLE) zamverb$(LIB_EXT)

.PHONY: clean install uninstall
