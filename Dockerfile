FROM ubuntu:trusty


RUN apt-get update && apt-get install -y gcc make vim  python


WORKDIR /usr/src/app/
#COPY * /usr/src/app/


