DOCKER_RUN = docker run --rm -v $(shell pwd):/src -u $(shell id -u):$(shell id -g)

EXPORT_METHODS = #-sEXPORTED_FUNCTIONS='["feval"]'
EMFLAGS = $(EXPORT_METHODS) -sFILESYSTEM=0 -sALLOW_MEMORY_GROWTH=1 -sSINGLE_FILE=1 # --closure 1 #-sASSERTIONS=2 
CXXFLAGS = -I/usr/include/eigen3/ -O3 -lembind -Wall -Werror -fno-exceptions $(EMFLAGS)

aaa.js: aaa.cpp
	$(DOCKER_RUN) -t aaa em++ $(CXXFLAGS) aaa.cpp -o aaa.js --post-js postamble.js
	mv aaa.js ../build/aaa.js
#mv aaa.wasm ../build/aaa.wasm

build: aaa.js

clean:
	rm -f ../build/aaa.js ../build/aaa.wasm

build-docker:
	docker build -t aaa .
