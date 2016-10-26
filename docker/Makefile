# Build Docker image
build: docker_build

# Image name can be overidden with env vars.
DOCKER_IMAGE ?= humburg/pyrates

# Use the version number as the release tag.
DOCKER_TAG = $(CODE_VERSION)

docker_build:
	# Build Docker image
	cp ../requirements.txt .
	docker build \
  --label org.label-schema.build-date=`date -u +"%Y-%m-%dT%H:%M:%SZ"` \
  --label org.label-schema.version=$(CODE_VERSION) \
  --label org.label-schema.vcs-url=`git config --get remote.origin.url` \
	-t $(DOCKER_IMAGE) .