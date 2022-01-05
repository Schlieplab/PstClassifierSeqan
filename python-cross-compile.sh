docker build --rm -t pst-classifier-seqan-python-cross-compile -f python-package/Dockerfile .

echo "Copying artifacts from docker container"
docker run --rm -v "$PWD/python-package/dist":/out/dist -t pst-classifier-seqan-python-cross-compile:latest
