docker build -t dockerfile . |& tee docker.txt

docker-machine stop
docker-machine start

docker rm $(docker ps -a -q)
docker rm $(docker ps -a -q)

docker run -e PASSWORD=pass -p 8787:8787 dockerfile
