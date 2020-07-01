FROM bucket.canfar.net/gem2caom2

RUN pip install matplotlib

ARG OMC_REPO=opencadc-metadata-curation

RUN cd /tmp && \
  git clone https://github.com/${OMC_REPO}/caom2pipe.git && \
  pip install ./caom2pipe

RUN cd /tmp && \
  git clone https://github.com/${OMC_REPO}/gem2caom2.git && \
  pip install ./gem2caom2

RUN git clone https://github.com/${OMC_REPO}/gemProc2caom2.git && \
  cp ./gemProc2caom2/scripts/config.yml / && \
  cp ./gemProc2caom2/scripts/docker-entrypoint.sh / && \
  pip install ./gemProc2caom2

ENTRYPOINT ["/docker-entrypoint.sh"]
