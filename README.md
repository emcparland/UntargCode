# Pipeline for pre-processing a multi-batch untargeted exometabolome experiment with XCMS on a HPC

*A big thank you to Krista Longnecker (WHOI) who laid the groundwork for this code and Elzbieta Lauzikaite (Imperial College London) who setup [a similar framework for pbs](https://github.com/lauzikaite/Imperial-HPC-R) that I built off*

**As is, you should be able to run this in your own compute space by installing the conda environment and altering the paths and inputs in the run files. But if you have a different experimental setup, you should also have a look at the R scripts**

## Install the conda environment via the yml file:
```conda env create --file untargmetab.yml```

This includes R version 3.6 plus XCMS3 and Autotuner, and jupyter notebook for later analyses. If you're not comfortable with conda or conda+R I recommend starting by reading  this [blog post by Sarah Hu](https://alexanderlabwhoi.github.io/post/anaconda-r-sarah/) and then use your friend google.

### Note about activating this conda environment on hpc with slurm:
- Remember that your sbatch command will create a new compute environment for each array so it doesn't know about your conda init
- Most scripts I've seen use the 'source activate myenv' command. However, when using this my .Rout file kept randomly throwing this error: "/vortexfs1/home/emcparland/.conda/envs/untargmetab/lib/R/bin/R: line 240: /vortexfs1/home/emcparland/.conda/envs/untargmetab/lib/R/etc/ldpaths: No such file or directory"
- This seems to be a problem for [others](https://github.com/conda-forge/r-base-feedstock/issues/67)
- One suggested solution was to use the newer command 'conda activate myenv'. To do this you need to source your conda.sh, you'll see this in the code as:

```CONDA_BASE=$(conda info --base)```
```source $CONDA_BASE/etc/profile.d/conda.sh```
```conda activate untargmetab```

-However, the script was still giving the same error randomly when I run larger arrays (though not as many as before?)

-One step further, seems like this is an issue of running the array and initializing the environment every time. For some reason the conda activate initializes the path every time and sometimes the path doesn't exist? I don't fully understand this yet but it seems to be an unresolved issue on git. A fix proposed by another git user and that seems to be working for me is to edit the activate-r-base.sh script in the environment:

```nano ~/.conda/envs/untargmetab/etc/conda/activate.d/activate-r-base.sh```

comment out the "R CMD javareconf" line to look like this: 

```#!/usr/bin/env sh```

```#R CMD javareconf > /dev/null 2>&1 || true ```

## Step 1: Create metadata
This is a quick R script to create a tab-delimited metadata file of all the sequence files (if you have multiple batches) and keep only the mzML files you want to peak pick and align (e.g. I remove the 9 conditioning pool samples here from each batch). *Make sure you have added a column named ionMode (pos or neg) and goodData (0 or 1).* It will also add an extra column to the metadata with the path of each mzml file that is useful for later.

```sbatch scripts_dir/run-metadata.slurm```

Check how many files you have 
```wc -l metadata.txt```

I have 502 and I will use this number in Step 3 to set the total number of array jobs that will be run.

## Step 2: peak picking and peak shape evaluation
Run the peak picking and peak shape on each file individually with an array job. I chose my parameters based on optimization of picking our 22 stable isotope labeled internal standards. This step is an 'embarassingly parallel' computation so I use a job array to quickly process hundreds of files. I run 40 jobs at a time and each jobs takes about 20 minutes each. I filter the peaks based on RMSE < 0.125 Then use peak cleaning functions to remove wide peaks (<40 s) and merge neighboring peaks. For 500 files, I am done with Step 3 in ~3 hours :clap: :grin: :clap:

```sbatch scripts_dir/run-xcms1.slurm```

Update status of jobs to your screen if you're interested (this is how I discovered the issue mentioned above of skipping files) ```watch -n 60 squeue -u emcparland```

## Step 3: combine picked peaks
To speed up peak picking, we performed peak picking as an array. Now combine into a single MS OnDisk object

```sbatch scripts_dir/run-xcms_combine.slurm```

## Step 4: perform retention time correction, grouping and fill peaks
This will use xcms to clean up peak picking with refineChromPeaks, then perform orbiwarp retention time correction, correspondence (peak grouping), and fill peaks. As I ran a pooled sample every five samples in these batches, I use the subset option for retention time alignment and peak grouping. At each stage a new RData object is saved in case something crashes in the middle or you want to look at the files while they are running. Finally it will output two csv files, one with all of the peaks ("aligned.csv") and the second with the feature count table ("picked.csv")

Note: For reference, when I was testing this code with ~100 samples, I could run this on one 'small' memory node of 185GB. However, my actual dataset of 500+ samples required being run on the 'bigmem' partition with 500GB of memory. The refinechrompeaks and fill peaks steps require loading the original raw files and therefore required the bigmem memory space (obiwarp and correspondence require much less memory).

```sbatch scripts_dir/run-xcms2.slurm```

## Step 5: Create an xset object 
CAMERA will require your data object to be in the 'old' XCMS format. This script will create this object for you. Note the fix-around for the error thrown by sample class naming. I had to use bigmem to make fillPeaks run. Make sure you edit the polarity mode.

```srun -p bigmem --time=04:00:00 --ntasks-per-node=1 --mem=500gb --pty bash```

```conda activate untargmetabR4```

```R```

```source("create_xset.R")```

## Step 6: Use CAMERA to create pseudospectra.
CAMERA also uses the xcmsSet object. You will need both positive and negative ionization mode output here.

```sbatch scripts_dir/run-camera.slurm```


## Stable isotope labeled internal standards.
All samples were spiked with a mix of the following internal standards.

| compound | label | final ng/ml | ion mode | exact mass |
| -------- | ----- | ----------- | -------- | ---------- |
| 3m2ob    | 13C5  | 150 | negative | 121.0641 |
| alpha-ketoglutaric acid | D6 | 2000 | negative | 152.0592 |
| betaine | D11 | 20 | positive | 128.1481 |
| citric acid | D4 | 250 | negative | 196.0521 |
| cysteine | D3 | 1000 | positive | 124.0385 |
| guanosine | D2 | 50 | both | 285.1043 |
| malic acid (DL) | D3 | 100 | negative | 137.0403 |
| phenyalanine | D8 | 50 | both | 173.1292 |
| proline | 13C5,15N | 25 | positive | 121.0771 |
| lysine | D4 | 750 | positive | 150.1306 |
| xanthine | 15N2 | 50 | both | 154.0275 |
| 4 aminobenzoic acid | D4 | 2000 | positive | 141.0728 |
| sn-glycerol-3-phosphate | 13C3 | 500 | both | 175.0238 |
| pyridoxine (vitamin B6) | 13C4 | 25 | both | 173.0873 |
| pantothenate | 13C3,15N | 25 | both | 223.1178 |
| leucine | D3 | 50 | both | 134.1134 |
| AMP | 15N5 | 100 | both | 352.0483 |
| methionine | D3 | 50 | both | 152.0698 |
| biotin | D2 | 50 | both | 246.1008 |
| succinic acid | D6 | 500 | negative | 124.0643 |
| cholic acid | D4 | 50 | negative | 412.3127 |
| indole-3-acetic acid | D7 | 1000 | negative | 182.1072 |
