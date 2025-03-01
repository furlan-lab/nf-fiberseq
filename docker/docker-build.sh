#!/bin/bash

cd ./docker/pacbio
docker buildx build --platform linux/amd64 -t tamdbho/pacbio:v1.3 -t tamdbho/pacbio:latest --push .

cd ./docker/fiberseq
docker buildx build --platform linux/amd64 -t tamdbho/fiberseq:v1.3 -t tamdbho/fiberseq:latest --push .

cd /Users/tho3/Desktop/docker/fiberseq-conda
docker buildx build --platform linux/amd64 -t tamdbho/fire:v1.0 -t tamdbho/fire:latest --push .
