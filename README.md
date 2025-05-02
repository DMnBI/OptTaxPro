## OptTaxPro: Optimized Taxonomic Profiling of microbiome using full-length 16S rRNA sequences

This repository includes the implementation and experimental data of 'OptTaxPro: Optimized Taxonomic Profiling of microbiome using full-length 16S rRNA sequences'. Please cite our paper if you use our pipeline. Fill free to report any issue for maintenance of our model.

## Citation
If you have used OptTaxPro in your research, please cite the following publication:

*In Review*


## 1. Setup
We strongly recommend you to use python virtual environment with [Anaconda](https://www.anaconda.com)/[Miniconda](https://docs.conda.io/en/latest/miniconda.html). The details of the environment we used are as follows:

* python 3.9.18
* cutadapt 5.0
* ete3 3.1.3
* numpy 1.26.4
* pandas 2.2.3
* vsearch 2.30.0
* scikit-learn 1.2.2 (only for building custom HSGs)
* scikit-learn-extra 0.3.0 (only for building custom HSGs)


### 1.1 Build environment
```
conda update conda (optional)
cd OptTaxPro/
conda env create -f environment.yml -n OptTaxPro
conda activate OptTaxPro
```

Please make sure all the required programs are successfully installed


### 1.2 Decompress databases

The source files and useful scripts are in this repository. The database files have been uploaded on **data** directory. Decompress database files before use

```
cd data/
gunzip *.gz
```

### 1.3 Allow executable permissions

For convenience, you need to allow executable permissions for all scripts

```
# On the top of the directory
chmod +x OptTaxPro/OptTaxPro
find . -name '*.py' -type f | xargs chmod +x
```


## 2. How to use OptTaxPro

OptTaxPro consists of **FOUR** main functions: **preprocess**, **cluster**, **classify**, **profile**  
There is a main wrapper script `OptTaxPro` in the `OptTaxPro` directory  
In addtion, there is an end-to-end function: **alltheway**

```
OptTaxPro {preprocess, cluster, classify, profile, alltheway} [options]
```
you can find details of required/optional parameters for each function with -h option.

```
OptTaxPro {preprocess, cluster, classify, profile, alltheway} -h
```

### 2.1 Preprocess
The input sequences are being filtered through adaptor removal, quality control, and singleton removal sequentially

```
export DATA=/path/to/your/data
export PRE=/path/to/filtered/reads
export NUM_CPU="the number of CPUs will be used"

# basic command
./OptTaxPro preprocess -i $DATA -o $PRE

# using $NUM_CPU threads
./OptTaxPro preprocess -i $DATA -o $PRE -p $NUM_CPU

# saving log file
./OptTaxPro preprocess -i $DATA -o $PRE -p $NUM_CPU --log preprocessing.log

# adjusting sequence length cutoffs
./OptTaxPro preprocess -i $DATA -o $PRE --min-len 1000 --max-len 1600

# adjusting sequence quality cutoff
./OptTaxPro preprocess -i $DATA -o $PRE --qc-cutoff 0.996

# allowing only exactly matched primers
./OptTaxPro preprocess -i $DATA -o $PRE -e 0.0

# adjusting singleton filtering cutoff
./OptTaxPro preprocess -i $DATA -o $PRE --single-cutoff 0.99

# allowing more threads for clustering
./OptTaxPro preprocess -i $DATA -o $PRE --single-threads 10
```

### 2.2 Clustering
OTUs are going to be built using VSEARCH

```
export OTU=/path/to/OTUs

# basic command
./OptTaxPro cluster -i $PRE -o $OTU --suffix .clean.fastq

# saving clustering log
./OptTaxPro cluster -i $PRE -o $OTU --suffix .clean.fastq --log cluster.log

# adjusting clustering threshold
./OptTaxPro cluster -i $PRE -o $OTU --suffix .clean.fastq --otu-cutoff 0.97

# allowig more threads for clustering
./OptTaxPro cluster -i $PRE -o $OTU --suffix .clean.fastq --otu-threads 10
```

### 2.3 Classify
Taxonomy assignment is being conducted by homology search

```
export CLASSIFY=/path/to/classification/results

# basic command
./OptTaxPro classify -i $OTU -o $CLASSIFY

# using $NUM_CPU threads
./OptTaxPro classify -i $OTU -o $CLASSIFY -p $NUM_CPU

# saving classification logs
./OptTaxPro classify -i $OTU -o $CLASSIFY -p $NUM_CPU --log classify.log

# using custom DB
./OptTaxPro classify -i $OTU -o $CLASSIFY -p $NUM_CPU --db {PATH_TO_CUSTOM_DB}

# adjusting intermediate search cutoffs
./OptTaxPro classify -i $OTU -o $CLASSIFY -p $NUM_CPU --search-cutoffs 0.97 0.94 0.85 0.75 0.7

# adjusting assignment ranks
./OptTaxPro classify -i $OTU -o $CLASSIFY -p $NUM_CPU --ranks species genus family

# adjusting assignment cutoffs
./OptTaxPro classify -i $OTU -o $CLASSIFY -p $NUM_CPU --ranks species genus family --assign-cutoffs 99 94 86

# appending scientific_name column into output table
./OptTaxPro classify -i $OTU -o $CLASSIFY -p $NUM_CPU --add-name

# expanding results for all queries
./OptTaxPro classify -i $OTU -o $CLASSIFY -p $NUM_CPU --u_dir $OTU --uc-suffix .uc
```

Columns of `taxonomy assignment` file:

|Columns | Description |
|:-------|:------------|
|seqid |Query sequence ID |
|rank |classified rank |
|taxid |NCBI taxonomy ID of besthit(s) |
|identity |identity of besthit(s) |
|assigned | **taxid** if identity > cutoff (rank) else **-1** |
|scientific_name | translated assigned taxid (unclassified for -1) |

**NOTE**: 'scientific_name' column will only be added when '--add-name' option is given

### 2.4 Profile
Estimating taxonomy profile based on the taxonomy assignment results

```
export PROFILE=/path/to/profiling/outputs

# basic command
./OptTaxPro profile -i $CLASSIFY -o $PROFILE

# using $NUM_CPU threads
./OptTaxPro profile -i $CLASSIFY -o $PROFILE -p $NUM_CPU

# saving classification logs
./OptTaxPro profile -i $CLASSIFY -o $PROFILE -p $NUM_CPU --log profile.log

# convert taxid into scientific name
./OptTaxPro profile -i $CLASSIFY -o $PROFILE -p $NUM_CPU --base-col scientific_name

# report additional ranks
./OptTaxPro profile -i $CLASSIFY -o $PROFILE -p $NUM_CPU --base-col scientific_name --profile-ranks species HSG genus family

# apply multiple filtering cutoffs 
./OptTaxPro profile -i $CLASSIFY -o $PROFILE -p $NUM_CPU --filtering-cutoffs 0 0.5 1.5

# adjusting filtering method
./OptTaxPro profile -i $CLASSIFY -o $PROFILE -p $NUM_CPU --filtering-pivot mean

# modify prefix of output tables
./OptTaxPro profile -i $CLASSIFY -o $PROFILE -p $NUM_CPU --output-prefix mysample
```

### 2.5 End-to-end RUN
Users can simply run all the processes using one command

```
export DATA=/path/to/your/data
export OUTPUT=/path/to/profiling/outputs
export NUM_CPU="the number of CPUs will be used"

# basic command
./OptTaxPro alltheway -i $DATA -o $OUTPUT --base-col scientific_name --output-prefix mysample

# using #NUM_CPU threads
./OptTaxPro alltheway -i $DATA -o $OUTPUT --base-col scientific_name --output-prefix mysample -p $NUM_CPU

# saving all logs
./OptTaxPro alltheway -i $DATA -o $OUTPUT --base-col scientific_name --output-prefix mysample -p $NUM_CPU --log alltheway.log
```

### 2.6 Using configuration file
Users can adjust and save all optional parameters using configuration file. An example of the configuration file is included in **data** directory. `data/config.cfg`

```
# preprocess
./OptTaxPro preprocess -c config.cfg

# cluster
./OptTaxPro cluster -c config.cfg

# classify
./OptTaxPro classify -c config.cfg

# profile
./OptTaxPro profile -c config.cfg

# alltheway
./OptTaxPro alltheway -c config.cfg
```

## 3. Build custom HSG

We provide a script for building custom **Homologous Species Group (HSG)**. The scripts are placed in the `scripts` directory.

### 3.1 build_HSG.py

This script conducts HSG building algorithm that takes 16S sequences as an input and an HSG table as output.

**NOTE1.** All sequences have to include NCBI accession number in the beginning of their header  
(This script automatically detects and splits by their genus according to it) 

**USAGE**

```
build_HSG.py \
    -i examples/build_HSGs/three_genera.fna \
    -a OptTaxPro/data/acc2taxid.txt \
    -o three_genera.HSG.csv
```

## 4. Run example data

This is an example of OptTaxPro workflow to classify the simulated datasets. Example data are given in the `examples/OptTaxPro` directory. This is the same simulated sequences that used in the original article.

### 4.1 RUN OptTaxPro

**NOTE** These simulated reads are all preprocessed ones. Thus, you **DO NOT** run preprocessing step for this (all reads are going to be ignored due to lack of primer sequences) 

```
# making temporary directory
mkdir test_run

# Building OTUs first
./OptTaxPro/OptTaxPro cluster \
	-i examples/OptTaxPro \
	-o test_run/OTUs \
	--suffix .pbsim.fasta \
	--log test_run.log
	
# Performing Taxonomic assignments
./OptTaxPro/OptTaxPro classify \
	-i test_run/OTUs \
	-o test_run/classify \
	--t_dir test_run \
	--log test_run.log \
	--search-cutoffs 0.97 \
	--ranks species HSG genus \
	--assign-cutoffs 97 97 94 \
	--remove-self \
	--u_dir test_run/OTUs \
	--add-name
	
# Making profile table
./OptTaxPro/OptTaxPro profile \
	-i test_run/classify \
	-o test_run/profile \
	--log test_run.log \
	--profile-ranks species HSG genus \
	--base-col scientific_name \
	--output-prefix simulated
```

### 4.2 RUN OptTaxPro for realworld data

**NOTE** Before running OptTaxPro, download samples from NCBI SRA under the accession number (PRJNA933120) and save them in `examples/PRJNA933120` directory.

```
# making temporary directory
mkdir test_realworld

# Run end-to-end process using config file
./OptTaxPro/OptTaxPro alltheway \
	-c OptTaxPro/data/config.cfg 
```

### 4.3 Build custom HSGs

```
# On the top of the directory

./scripts/build_HSG.py \
	-s examples/build_HSGs/three_genera.fna \
	-a OptTaxPro/data/acc2taxid.txt \
	-o HSG_three_genera.csv
```

