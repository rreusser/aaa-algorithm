FROM emscripten/emsdk

RUN apt-get update
RUN apt-get install -y libeigen3-dev

WORKDIR /src
RUN mkdir -p /src
VOLUME ["/src"]
CMD ["/bin/bash"]