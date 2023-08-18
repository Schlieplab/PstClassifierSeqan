FROM python:3.10-bullseye as python310

COPY src /PstClassifierSeqan/src
COPY include /PstClassifierSeqan/include
COPY seqan3 /PstClassifierSeqan/seqan3
COPY robin-hood-hashing /PstClassifierSeqan/robin-hood-hashing
COPY HighFive /PstClassifierSeqan/HighFive
COPY eigen /PstClassifierSeqan/eigen
COPY vlmc-from-kmers /PstClassifierSeqan/vlmc-from-kmers
COPY indicators /PstClassifierSeqan/indicators
COPY python-package/ /PstClassifierSeqan/python-package/
COPY python-package/libvlmc /PstClassifierSeqan/python-package/libvlmc
COPY python-package/setup.py /PstClassifierSeqan/python-package/setup.py
COPY python-package/CMakeLists.txt /PstClassifierSeqan/python-package/CMakeLists.txt
COPY python-package/dev-requirements.txt /PstClassifierSeqan/python-package/dev-requirements.txt

RUN apt-get update -y && apt-get install -y cmake ninja-build libhdf5-dev libeigen3-dev libtbb-dev libgomp1 doxygen


WORKDIR /PstClassifierSeqan/python-package/

RUN python3.10 -m pip install -r dev-requirements.txt
RUN python3.10 setup.py bdist_wheel

FROM python:3.9-bullseye as python39

COPY src /PstClassifierSeqan/src
COPY include /PstClassifierSeqan/include
COPY seqan3 /PstClassifierSeqan/seqan3
COPY robin-hood-hashing /PstClassifierSeqan/robin-hood-hashing
COPY HighFive /PstClassifierSeqan/HighFive
COPY eigen /PstClassifierSeqan/eigen
COPY vlmc-from-kmers /PstClassifierSeqan/vlmc-from-kmers
COPY indicators /PstClassifierSeqan/indicators
COPY python-package/ /PstClassifierSeqan/python-package/
COPY python-package/libvlmc /PstClassifierSeqan/python-package/libvlmc
COPY python-package/setup.py /PstClassifierSeqan/python-package/setup.py
COPY python-package/CMakeLists.txt /PstClassifierSeqan/python-package/CMakeLists.txt
COPY python-package/dev-requirements.txt /PstClassifierSeqan/python-package/dev-requirements.txt

RUN apt-get update -y && apt-get install -y cmake ninja-build libhdf5-dev libeigen3-dev libtbb-dev libgomp1 doxygen


WORKDIR /PstClassifierSeqan/python-package/

RUN python3.9 -m pip install -r dev-requirements.txt
RUN python3.9 setup.py bdist_wheel

FROM python:3.8-bullseye as python38

COPY src /PstClassifierSeqan/src
COPY seqan3 /PstClassifierSeqan/seqan3
COPY robin-hood-hashing /PstClassifierSeqan/robin-hood-hashing
COPY HighFive /PstClassifierSeqan/HighFive
COPY eigen /PstClassifierSeqan/eigen
COPY vlmc-from-kmers /PstClassifierSeqan/vlmc-from-kmers
COPY indicators /PstClassifierSeqan/indicators
COPY python-package/ /PstClassifierSeqan/python-package/
COPY python-package/libvlmc /PstClassifierSeqan/python-package/libvlmc
COPY python-package/setup.py /PstClassifierSeqan/python-package/setup.py
COPY python-package/CMakeLists.txt /PstClassifierSeqan/python-package/CMakeLists.txt
COPY python-package/dev-requirements.txt /PstClassifierSeqan/python-package/dev-requirements.txt

RUN apt-get update -y && apt-get install -y cmake ninja-build libhdf5-dev libeigen3-dev libtbb-dev libgomp1 doxygen


WORKDIR /PstClassifierSeqan/python-package/

RUN python3.8 -m pip install -r dev-requirements.txt
RUN python3.8 setup.py bdist_wheel

FROM python:3.7-bullseye as python37

COPY src /PstClassifierSeqan/src
COPY seqan3 /PstClassifierSeqan/seqan3
COPY robin-hood-hashing /PstClassifierSeqan/robin-hood-hashing
COPY HighFive /PstClassifierSeqan/HighFive
COPY eigen /PstClassifierSeqan/eigen
COPY vlmc-from-kmers /PstClassifierSeqan/vlmc-from-kmers
COPY indicators /PstClassifierSeqan/indicators
COPY python-package/libvlmc /PstClassifierSeqan/python-package/libvlmc
COPY python-package/setup.py /PstClassifierSeqan/python-package/setup.py
COPY python-package/CMakeLists.txt /PstClassifierSeqan/python-package/CMakeLists.txt
COPY python-package/dev-requirements.txt /PstClassifierSeqan/python-package/dev-requirements.txt

RUN apt-get update -y && apt-get install -y cmake ninja-build libhdf5-dev libeigen3-dev  libtbb-dev libgomp1 doxygen

WORKDIR /PstClassifierSeqan/python-package/

RUN python3.7 -m pip install -r dev-requirements.txt
RUN python3.7 setup.py bdist_wheel


FROM python:3.6-bullseye as python36

COPY src /PstClassifierSeqan/src
COPY seqan3 /PstClassifierSeqan/seqan3
COPY robin-hood-hashing /PstClassifierSeqan/robin-hood-hashing
COPY HighFive /PstClassifierSeqan/HighFive
COPY eigen /PstClassifierSeqan/eigen
COPY vlmc-from-kmers /PstClassifierSeqan/vlmc-from-kmers
COPY indicators /PstClassifierSeqan/indicators
COPY python-package/libvlmc /PstClassifierSeqan/python-package/libvlmc
COPY python-package/setup.py /PstClassifierSeqan/python-package/setup.py
COPY python-package/CMakeLists.txt /PstClassifierSeqan/python-package/CMakeLists.txt
COPY python-package/dev-requirements.txt /PstClassifierSeqan/python-package/dev-requirements.txt

RUN apt-get update -y && apt-get install -y cmake ninja-build libhdf5-dev libeigen3-dev libtbb-dev libgomp1 doxygen

WORKDIR /PstClassifierSeqan/python-package/

RUN python3.6 -m pip install -r dev-requirements.txt
RUN python3.6 setup.py bdist_wheel

FROM python:3.10-bullseye

WORKDIR /artifacts/dist/

COPY --from=python310 /PstClassifierSeqan/python-package/dist/libvlmc-* /artifacts/dist/
COPY --from=python39 /PstClassifierSeqan/python-package/dist/libvlmc-* /artifacts/dist/
COPY --from=python38 /PstClassifierSeqan/python-package/dist/libvlmc-0.2-cp38-cp38-linux_x86_64.whl /artifacts/dist/
COPY --from=python37 /PstClassifierSeqan/python-package/dist/libvlmc-0.2-cp37-cp37m-linux_x86_64.whl /artifacts/dist/
COPY --from=python36 /PstClassifierSeqan/python-package/dist/libvlmc-0.2-cp36-cp36m-linux_x86_64.whl /artifacts/dist/

CMD cp /artifacts/dist/* /out/dist
