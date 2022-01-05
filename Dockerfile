FROM buildpack-deps:18.04

RUN apt-get update -y && \
    apt-get install -y \
    unzip \
    cmake \
    libhdf5-dev \
    libeigen3-dev \
    libboost-dev \
    libboost-system-dev \
    libboost-serialization-dev

WORKDIR /PstClassifierSeqan

COPY vlmc-from-kmers vlmc-from-kmers
COPY indicators indicators
COPY seqan3 seqan3
COPY eigen eigen
COPY robin-hood-hashing robin-hood-hashing
COPY HighFive HighFive

COPY src src
COPY tests tests
COPY CMakeLists.txt CMakeLists.txt

WORKDIR /PstClassifierSeqan/build

RUN cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make pst-classifier pst-batch-training pst-score-sequences bic

WORKDIR /PstClassifierSeqan/build/bin

RUN cp ../src/pst-classifier ../src/pst-batch-training ../src/pst-score-sequences ../src/bic /PstClassifierSeqan/build/bin

ENV PATH /PstClassifierSeqan/build/bin:$PATH

ENTRYPOINT ["bic"]
