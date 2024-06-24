# Imaging Clusters of Pediatric Low-Grade Glioma are Associated with Distinct Molecular Characteristics

This repository includes codes relevant to the manuscript "Imaging Clusters of Pediatric Low-Grade Glioma are Associated with Distinct Molecular Characteristics" by 

"Anahita Fathi Kazerooni*, Adam Kraya*, Komal S. Rathi, Meen Chul Kim, Varun Kesherwani, Ryan Corbett,  Arastoo Vossough, Nastaran Khalili, Deep Gandhi, Neda Khalili, Ariana M. Familiar, Run Jin, Xiaoyan Huang, Yuankun Zhu, Alex Sickler, Matthew R. Lueder, Saksham Phul, Mateusz Koptyra, Yuanquan Song, Phillip B. Storm, Christos Davatzikos, Jeffrey B. Ware, Jessica Foster, Sabine Mueller, Jo Lynne Rokita, Michael J. Fisher, Adam C. Resnick, Ali Nabavizadeh" 

*Equally-contributing authors

## Imaging Clustering
### Software Requirements
- CaPTk, v1.8.1 (https://cbica.github.io/CaPTk/)
- Python3

### Hardware Used for this Study
- CUBIC (HPC Cluster) (https://www.med.upenn.edu/cbica/cubic.html)
- AWS/EC2 for batch image pre-processing, segmentation, radiomic feature extraction


### MRI Pre-processing and Tumor Segmentation:
- Required MRI sequences: T1, T1CE, T2, FLAIR (ADC optional)
- Pre-processing using BraTS Pre-processing Pipeline, details explained in: https://cbica.github.io/CaPTk/preprocessing_brats.html
- Tumor Segmentation, all details provided in: https://github.com/d3b-center/peds-brain-auto-seg-public
- Skull-stripping to generate a brain mask: https://github.com/d3b-center/peds-brain-auto-skull-strip
- Image normalization:
   - run_rescale.py
- Whole tumor generation:
   - run_wtmask.py
- Radiomic feature extraction, using CaPTk v1.8.1: https://cbica.github.io/CaPTk/ht_FeatureExtraction.html;
  - parameter file for radiomic feature extraction: radiomic_feature_params_20230725.csv
  - sample batch file: SampleBatchFile.csv

### Clustering:

## Transcriptomic Analysis
page

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
