#!/bin/bash
# Build the container image and drop into a shell with the repo mounted at /app.
set -e
cd "$(dirname "$0")/.."
docker build -t atmosphere-model:latest .
docker run --rm -it --name atmosphere-model -v "$PWD:/app" atmosphere-model:latest
