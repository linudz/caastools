# CAAStools Dockerization 1.0 - Software documentation.

## CAASTools: Docker barebones containerization WIP

### CAASTools Barebones Container WIP
- To build the Docker Image:
`docker build -t <image_name> .`  # Estimated build time around 4m30s (due to RERConverge)

- Create a container: 
`docker run -it --name caastools-barebones caastools-barebones`

- CT is already added to $PATH and so can be used through docker exec, or be introduced as parameter in workflow tools (such as Nextflow)

#### RERConverge docker instance
Run from the official docker container. RStudio enabled, can also be used through CLI and can be piped through workflow tools (I guess)
`docker run -it -p 8787:8787 -e PASSWORD=<password> --name <container_name> wem26/rerconverge`

## PLEASE NOTE THAT THE TOOL HAS NOT BEEN TRIED IN CONTAINER MODE AND FURTHER MODIFICATIONS ARE MOST PROBABLY NEEDED, THIS IS A WIP
