# Dockerfile to build container for unit testing

FROM openjdk:10

RUN apt-get update && apt-get install -y git ant

ENV DISPLAY :0
ENV USER root

WORKDIR /root

ADD . ./

ENTRYPOINT ant test
