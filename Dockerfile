FROM buildpack-deps:18.04

RUN apt-get update -y && \
    apt-get install -y \
    unzip \
    cmake

WORKDIR /PstClassifierSeqan

COPY seqan3 seqan3
COPY src src
COPY tests tests

WORKDIR /PstClassifierSeqan/build

RUN cmake -DCMAKE_BUILD_TYPE=Release ../src
RUN make

ENV PATH /PstClassifierSeqan/build:$PATH

WORKDIR /PstClassifierSeqan/tests_build

RUN cmake -DCMAKE_BUILD_TYPE=Release ../tests
RUN make

CMD ./probabilistic_suffix_tree_tests
