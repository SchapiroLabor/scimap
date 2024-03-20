FROM python:3.8

LABEL maintainer="Miguel Ibarra; Chiara Schiller, Jose Nimo"

ENV MPLCONFIGDIR=/tmp
ENV NUMBA_CACHE_DIR=/tmp
RUN chmod 777 /tmp

RUN apt-get update -qq && apt-get install -y \
    build-essential \
    ffmpeg \
    libsm6 \
    libxext6

RUN pip install --upgrade pip
RUN pip install --no-cache-dir scimap --upgrade

WORKDIR /scimap

COPY . .
