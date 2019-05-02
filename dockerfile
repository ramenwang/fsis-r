FROM scratch
USER root
RUN chmod +x /tmp
ARG GITHUB_TOKEN = ""
ENV GITHUB_PAT ${GITHUB_TOKEN}
RUN R -e 'install.packages(c("parallelDist", "spatstat", "FNN", "raster"), repos = "http://cran.wustl.edu", dependencies = T)'
COPY main /main
COPY shared /shared
WORKDIR /main
