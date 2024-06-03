# D3b TIRU Radiomic Analysis Repository

This repository includes codes relevant to each radiomic analysis project.  

### To reproduce the code in the `braf-fusions` module of this repository

1. Clone the repository:
```
git clone git@github.com:d3b-center/TIRU_radiomic_analysis.git
```

2. Pull Docker container:
```
docker pull pgc-images.sbgenomics.com/corbettr/tiru-radiomics:latest
```

3. Start the Docker container; from the `TIRU_radiomic_analysis` folder, run:
```
docker run --name <CONTAINER_NAME> -d -e PASSWORD=pass -p 8787:8787 -v $PWD:/home/rstudio/TIRU_radiomic_analysis pgc-images.sbgenomics.com/corbettr/tiru-radiomics:latest
```

NOTE: if using a Macbook with M1 chip, add the option `--platform linux/amd64` to above code

4. Execute the shell within the docker image; from the `TIRU_radiomic_analysis` folder, run: 
```
docker exec -ti <CONTAINER_NAME> bash
```
