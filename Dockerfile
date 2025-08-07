FROM uwgac/primed-file-checks:0.7.0

RUN cd /usr/local && \
    git clone https://github.com/UW-GAC/primed-mvp-gsr.git
