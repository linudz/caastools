# CAAStools 1.0 - Software documentation.

## Docker minimal containerization

### Docker specific stuff
- To build the Docker Image:
`docker build -t <image_name> .`  # Estimated build time around 4m30s (due to RERConverge)
##### Probably better to just set to download the R dep if needed already through R console (thus through the script, it would be nice to have an R-devtools ready platform for it though)
- To run create a container: 
`docker run -it --name <container_name> <image_name>` # Esto tengo q ver que pasa exactamente al correr el container
Piensa sobre añadir puertecito 
