FROM python:3.8

LABEL maintainer="Miguel Ibarra; Chiara Schiller"

RUN apt-get update -qq && apt-get install -y \
    build-essential \
    ffmpeg \
    libsm6 \
    libxext6

RUN pip install --no-cache-dir scimap --upgrade

WORKDIR /scimap

COPY . .
