# Variables
DOCKERFILE_NAME := Dockerfile_dsOmics
IMAGE_NAME := brgelab/rock-omics
DSOMICS_VERSION := 1.0.19
R_VERSION := 4.4

# Build and push the image to Docker Hub
all: build-all push-all

build-all:
	docker build -f $(DOCKERFILE_NAME) --no-cache=true --pull -t="$(IMAGE_NAME):$(DSOMICS_VERSION)-R$(R_VERSION)" . && \
	docker tag $(IMAGE_NAME):$(DSOMICS_VERSION)-R$(R_VERSION) $(IMAGE_NAME):$(DSOMICS_VERSION) && \
	docker tag $(IMAGE_NAME):$(DSOMICS_VERSION)-R$(R_VERSION) $(IMAGE_NAME):latest

push-all:
	docker image push $(IMAGE_NAME):$(DSOMICS_VERSION)-R$(R_VERSION) && \
	docker image push $(IMAGE_NAME):$(DSOMICS_VERSION) && \
	docker image push $(IMAGE_NAME):latest

build:
	docker build --pull -f $(DOCKERFILE_NAME) -t $(IMAGE_NAME):$(DSOMICS_VERSION)-R$(R_VERSION) .
