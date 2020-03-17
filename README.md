# gbs
GBP pipeline for single-end or paired-end short reads (Illumina or IonTorrent).


## Installation
1. Create a virtual environment using conda
2. Activate environment
3. Clone repository
4. Install requirements
5. Test installation

```
conda create -n gbs python=3.6
conda activate gbs
git clone https://github.com/duceppemo/gbs
cd gbs
conda install --file requirement.txt
python3 gbs.py -h
```
## Usage
```
usage: gbs.py [-h] -r /reference_genome.fasta -i /input_folder/ -o
              /output_folder/ -a /left_adapters.fasta [-t 48] [-m 48] [-pe]
              [-se]

Reference-based genotyping-by-sequencing (GBS) pipeline

optional arguments:
  -h, --help            show this help message and exit
  -r /reference_genome.fasta, --reference /reference_genome.fasta
                        Reference genome for read mapping
  -i /input_folder/, --input /input_folder/
                        Folder that contains the fastq files
  -o /output_folder/, --output /output_folder/
                        Folder to hold the result files
  -a /left_adapters.fasta, --adapters /left_adapters.fasta
                        Fasta file with sequence of adapters to trim at the
                        beginning of reads
  -t 48, --threads 48   Number of CPU. Default is maximum CPU available(48)
  -m 48, --memory 48    Memory in GB. Default is 85% of total memory (459)
  -pe, --paired-end     Input fastq files are paired-end
  -se, --single-end     Input fastq files are single-end
  ```
  
