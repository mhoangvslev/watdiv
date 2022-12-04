FROM ubuntu:22.04

RUN apt update && apt install -y --no-install-recommends \
        curl \
		cmake \
		libboost-all-dev \
        wamerican \
        build-essential \
        python3 \
	&& rm -rf /var/lib/apt/lists/*

COPY . /opt/watdiv/
WORKDIR /opt/watdiv/

RUN make rebuild && make install
#ENTRYPOINT [ "python3 -m http.server" ]