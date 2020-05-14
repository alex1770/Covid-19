all:: getdata
deb:: 

CFLAGS:= -std=c++11 -Wall
ifneq ($(filter deb,$(MAKECMDGOALS)),)
  CFLAGS:= $(CFLAGS) -g
else
  CFLAGS:= $(CFLAGS) -O3
endif

.PHONY: getdata
getdata:
	./downloadscript.sh
