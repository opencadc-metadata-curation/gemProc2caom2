FROM bucket.canfar.net/gem2caom2

RUN pip install matplotlib

ARG OPENCADC_BRANCH=master
ARG OPENCADC_REPO=opencadc
ARG PIPE_BRANCH=master
ARG PIPE_REPO=opencadc

RUN pip install git+https://github.com/${OPENCADC_REPO}/caom2pipe@${OPENCADC_BRANCH}#egg=caom2pipe

RUN pip install git+https://github.com/${PIPE_REPO}/gemProc2caom2@${PIPE_BRANCH}#egg=gemProc2caom2

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
