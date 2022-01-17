FROM bucket.canfar.net/gem2caom2

RUN pip install matplotlib

ARG PIPE_BRANCH=master
ARG PIPE_REPO=opencadc

RUN pip install git+https://github.com/${PIPE_REPO}/gemProc2caom2@${PIPE_BRANCH}#egg=gemProc2caom2

RUN useradd --create-home --shell /bin/bash cadcops
USER cadcops

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
