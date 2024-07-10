#This is the working docker file for MacOS m3 that worked with my application
FROM python:3

RUN apt-get update \
    && apt-get install -y binutils libproj-dev gdal-bin python3-gdal python3-pip python3-fiona

RUN apt-get update && apt-get install -y gdal-bin \
     apt-utils libgdal-dev

RUN apt-get update && \
    apt-get install -y gcc python3-dev musl-dev \
    gdal-bin \
    python3-gdal \
    libgeos-dev \
    libffi-dev \
    proj-bin \
    libopenblas-dev \ 
    python3-scipy python3-numpy python3-pandas \
    netcat-traditional \
    nano 
   
RUN apt-get update && apt-get install -y gdal-bin \
     apt-utils libgdal-dev

RUN pip install pygdal=="`gdal-config --version`.*"

RUN pip install -I GDAL=="`gdal-config --version`.*"

COPY requirements.txt requirements.txt

# RUN pip install -r requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

RUN pip install pystac_client

RUN pip install odc-stac

RUN pip install rioxarray

WORKDIR /app

COPY . .

EXPOSE 9001

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "9001"]

@##docker compose file is also here
services:
  app:
    build:
      context: .
      dockerfile: Dockerfile
    ports:
      - "9001:9001"
    volumes:
      - .:/app
    environment:
      - PYTHONUNBUFFERED=1
    restart: always

