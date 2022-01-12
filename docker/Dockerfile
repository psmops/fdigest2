FROM ubuntu:20.04
MAINTAINER schastel

RUN apt-get update && apt-get upgrade -y && apt-get install -y git make gcc wget
RUN mkdir -p /digest2 && cd /digest2 && git clone https://github.com/psmops/fdigest2 . && make

#
# docker build --tag=docker.io/schastel/psmops-digest2:0.19.2 .
# docker push docker.io/schastel/psmops-digest2:0.19.2
#
# docker run --rm -it docker.io/schastel/psmops-digest2:0.19.2 bash
# 