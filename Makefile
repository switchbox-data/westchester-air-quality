# End-to-end pipeline runner.
#
# Quickstart:
#   make data       # downloads the Aclima input data (see README)
#   export CENSUS_API_KEY=<your-key>
#   make            # builds image, runs 01a/01b/02/03/04
#
# Override knobs from the command line, e.g.:
#   make NMF_FACTORS=2,3,4,5,6   # full NMF sweep + IM/IS plot (~2h sequential)
#   make NMF_PARALLEL=4          # try parallel NMF; sometimes crashes in container
#   make 03                      # re-run a single stage

IMAGE        := westchester-air-quality:latest
NMF_FACTORS  ?= 3
NMF_PARALLEL ?= 1

# CENSUS_API_KEY must be exported in the caller's environment. Bail early
# with a clear message if it isn't.
ifndef CENSUS_API_KEY
$(warning CENSUS_API_KEY is not set — 02/03/04 will fail at the get_acs call)
endif

DOCKER_RUN = docker run --rm \
  -e CENSUS_API_KEY \
  -e NMF_FACTORS \
  -e NMF_PARALLEL \
  -v "$$PWD/Data:/project/Data" \
  -v "$$PWD/Results:/project/Results" \
  -v "$$PWD/Figures:/project/Figures" \
  $(IMAGE)

.PHONY: all data image 01 02 03 04 clean

all: image 01 02 03 04

data:
	scripts/fetch_data.sh

image:
	docker build -t $(IMAGE) .

# Script 01 is split per pollutant so each Rscript invocation peaks
# around 5–6 GB instead of 10–12 GB. See Code/01_Exposure_assessment.R.
01:
	$(DOCKER_RUN) bash -c 'POLLUTANTS=pm_2.5 Rscript Code/01_Exposure_assessment.R'
	$(DOCKER_RUN) bash -c 'POLLUTANTS=no2    Rscript Code/01_Exposure_assessment.R'

02:
	$(DOCKER_RUN) Rscript Code/02_Health_risk.R

03:
	$(DOCKER_RUN) Rscript Code/03_Source_apportionment.R

04:
	$(DOCKER_RUN) Rscript Code/04_Mitigation_strategies_NEW.R

clean:
	rm -rf Results/* Figures/*
