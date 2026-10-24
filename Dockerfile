# PanG2G container image: the g2g command plus the demo scripts and the test case.
#
#   docker build -t g2g .
#   docker run --rm g2g g2g --help
#   docker run --rm -v "$PWD":/data g2g g2g A.gfa B.gfa --ref GRCh38 --vcf delta.vcf
#   docker run --rm g2g /opt/g2g/container/selftest.sh        # the README's test case
#
# The engine is compiled in a build stage and copied into a slim runtime
# image (Ubuntu 24.04 with python3, numpy, matplotlib and pysam).  The binary
# targets x86-64-v3 (any AVX2 CPU) so the same image runs on any x86-64 node;
# add --build-arg G2G_MARCH=native to tune it to the build host instead.

FROM ubuntu:24.04 AS build
ARG G2G_MARCH=x86-64-v3
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
 && apt-get install -y --no-install-recommends g++ cmake make ca-certificates \
 && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/g2g
COPY . .
RUN rm -rf build \
 && cmake -B build/release -DCMAKE_BUILD_TYPE=Release -DG2G_BUILD_TESTS=OFF -DG2G_MARCH=${G2G_MARCH} \
 && cmake --build build/release -j"$(nproc)" --target g2g bench_matching \
 && build/release/g2g --help > /dev/null

FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        python3 python3-numpy python3-matplotlib python3-pysam \
        curl ca-certificates bzip2 \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /opt/g2g /opt/g2g
ENV PATH=/opt/g2g/build/release:${PATH} \
    G2G=/opt/g2g/build/release/g2g \
    MPLBACKEND=Agg
WORKDIR /data
CMD ["g2g", "--help"]
