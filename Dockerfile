FROM rocker/r-ver:4.5.1
RUN R -q -e "install.packages(c('Matrix','jsonlite','remotes'), repos='https://cloud.r-project.org')"
COPY . /opt/HeteroSpaSim
RUN R CMD INSTALL /opt/HeteroSpaSim
CMD ["R"]
