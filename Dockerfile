# rocker/geospatial only publishes multi-arch (incl. arm64) on :latest;
# the dated/versioned tags are amd64-only.
FROM rocker/geospatial:latest

# Install CRAN packages not already present in rocker/geospatial.
# rocker/geospatial already provides: sf, dplyr, ggplot2, lubridate, tidyr,
# reshape2, RColorBrewer, splines (base), tidycensus, tigris.
RUN install2.r --error --skipinstalled --ncpus -1 \
    DescTools \
    clipr \
    maptiles \
    ggnewscale \
    patchwork \
    tigris \
    tidycensus \
    showtext \
    sysfonts \
    svglite

# NMF depends on Bioconductor's Biobase, which is not on CRAN.
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org'); \
          BiocManager::install('Biobase', ask = FALSE, update = FALSE); \
          install.packages('NMF', repos='https://cloud.r-project.org')"

WORKDIR /project

# Create folders the R scripts write to. Data/ is expected to be bind-mounted.
RUN mkdir -p /project/Data /project/Results /project/Figures

COPY Code /project/Code

# Census ACS calls (tidycensus) need an API key. Pass via `docker run -e
# CENSUS_API_KEY=...` — get one at https://api.census.gov/data/key_signup.html.
CMD ["R", "--no-save"]
