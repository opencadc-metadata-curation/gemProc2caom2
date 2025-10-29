FROM opencadc/pandas:3.10-slim as builder

RUN apt-get update --no-install-recommends && \
    apt-get dist-upgrade -y && \
    apt-get install -y build-essential git && \
    rm -rf /var/lib/apt/lists/ /tmp/* /var/tmp/*

ARG CAOM2_BRANCH=master
ARG CAOM2_REPO=opencadc
ARG GEM_BRANCH=master
ARG GEM_REPO=opencadc-metadata-curation
ARG PIPE_BRANCH=master
ARG PIPE_REPO=opencadc-metadata-curation

RUN git clone https://github.com/${CAOM2_REPO}/caom2tools.git && \
    cd caom2tools && \
    git checkout ${CAOM2_BRANCH} && \
    pip install ./caom2 && \
    pip install ./caom2utils && \
    cd ..

RUN pip install git+https://github.com/${CAOM2_REPO}/caom2pipe@${CAOM2_BRANCH}#egg=caom2pipe

RUN pip install git+https://github.com/${GEM_REPO}/gem2caom2@${GEM_BRANCH}#egg=gem2caom2

RUN pip install git+https://github.com/${PIPE_REPO}/gemProc2caom2@${PIPE_BRANCH}#egg=gemProc2caom2

FROM python:3.10-slim
WORKDIR /usr/src/app

COPY --from=builder /usr/local/lib/python3.10/site-packages/ /usr/local/lib/python3.10/site-packages/
COPY --from=builder /usr/local/bin/* /usr/local/bin/
COPY --from=builder /usr/local/.config/* /usr/local/.config/

COPY --from=builder /etc/magic /etc/magic
COPY --from=builder /etc/magic.mime /etc/magic.mime
COPY --from=builder /usr/lib/x86_64-linux-gnu/libmagic* /usr/lib/x86_64-linux-gnu/
COPY --from=builder /usr/lib/file/magic.mgc /usr/lib/file/
COPY --from=builder /usr/share/misc/magic /usr/share/misc/magic
COPY --from=builder /usr/share/misc/magic.mgc /usr/share/misc/magic.mgc
COPY --from=builder /usr/share/file/magic.mgc /usr/share/file/magic.mgc

RUN useradd --create-home --shell /bin/bash cadcops
USER cadcops

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
