# OrthoTrees: Getting orthologs of the *Podospora* complex

A general pipeline to get orthologous genes from the *Podospora* complex. Due to the unpredictable number of ortholog groups, the pipeline relies on a [Snakemake checkpoint function](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html?fbclid=IwAR1v29DPDpWqve6yRlnc5vob2uIsxCfZt-NSjxfTtbaOZa4TFRuuqn8VbEk#data-dependent-conditional-execution).

## The data

The pipeline expects all the genome assemblies and annotations to be in a single folder, whose path is provided in the configuration file. These files have a name in the format `{strain}.nice.fa` and `{strain}.nice-XXX.gff3`, where `strain` stands for a specific strain code, and `XXX` is the version of the annotation. 

For reasons, I also keep the reference genomes in a different folder, so the configuration file expects the paths for those files. There are two references, Podan2 and PODCO. Each has its own gff3 file.

- Podan2: The reference genome of *P. anserina* (strain S), originally published by [Espagne et al. (2008)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2441463/), and improved in the Joint Genome Institute MycoCosm [website](https://genome.jgi.doe.gov/programs/fungi/index.jsf) under the name Podan2. There is a Podan3, but as far as I can tell is the same assembly. It's available also [here](https://github.com/johannessonlab/SpokPaper/blob/master/Fig4_S3_Backcrosses/extras/Podan2_AssemblyScaffoldsmt.fa).

- PODCO: The reference genome of *P. comata* is that of the strain T (also known as T<sub>D</sub> in [Vogan et al. 2019](https://elifesciences.org/articles/46454)) was published by [Silar et al. 2018](https://link.springer.com/article/10.1007/s00438-018-1497-3). It's deposited in the European Nucleotide Archive, [GCA_900290415.1](https://www.ebi.ac.uk/ena/data/view/GCA_900290415.1). Unfortunately the names are very ugly, so I modified them to be more like `Chromosome_1`, like such:

    $ cat GCA_900290415.1_version1_genomic.fna | sed 's;\(>[A-Z0-9\.]*\)\s\([a-zA-Z ,]*\): ;>Chromosome_;' > PODCO_genomic.fa
    $ sed -i 's;Chromosome_mitochondrion;Mitochondrion;' PODCO_genomic.fa

And for the gff file (less elegantly):

    $ sed -i 's/LR026964.1/Chromosome_1/' PODCO_genomic.gff3
    $ sed -i 's/LR026965.1/Chromosome_2/' PODCO_genomic.gff3
    $ sed -i 's/LR026966.1/Chromosome_3/' PODCO_genomic.gff3
    $ sed -i 's/LR026967.1/Chromosome_4/' PODCO_genomic.gff3
    $ sed -i 's/LR026968.1/Chromosome_5/' PODCO_genomic.gff3
    $ sed -i 's/LR026969.1/Chromosome_6/' PODCO_genomic.gff3
    $ sed -i 's/LR026970.1/Chromosome_7/' PODCO_genomic.gff3
    $ sed -i 's/LR026971.1/Mitochondrion/' PODCO_genomic.gff3

All reference files are available [here](https://github.com/johannessonlab/SpokPaper/tree/master/Fig1_3SppNetwork/references).

The pipeline will make symlinks of the data files in the working directory. 

## Extra scripts

The pipeline requires four external scripts, all available in [my GitHub](https://github.com/SLAment/Genomics). You have to provide their paths in the configuration file. They are:

- [gffutils2fasta.py](https://github.com/SLAment/Genomics/blob/master/GenomeAnnotation/gffutils2fasta.py)
- [orthogrs_parser.py](https://github.com/SLAment/Genomics/blob/master/Phylogenetics/orthogrs_parser.py)
- [query2hitseq.py](https://github.com/SLAment/Genomics/blob/master/BLAST/query2hitseq.py)
- [fastaconcat.py](https://github.com/SLAment/Genomics/blob/master/FastaManipulation/fastaconcat.py)

## Building the environment

First, I can start by updating conda.

    $ conda update -n base conda

Now, to create the environment.

    $ conda create -n LorePhylogenetics -c bioconda

**IMPORTANT!!** activate the environment before installing stuff! 
    
    $ conda activate LorePhylogenetics
    $ conda install -c bioconda snakemake-minimal=5.4.4
    $ conda install -c bioconda biopython=1.72=py37h04863e7_0
    $ conda install -c bioconda gffutils=0.9=py_1
    $ conda install -c bioconda mafft=7.407=1
    $ conda install -c bioconda iqtree=1.6.8

Unfortunately snakemake runs in python3 and OrthoFinder requires python 2. So I try to run the OrthoFinder rule with it's own conda environment.

    $ cat envs/orthofinder.yaml
    channels:
      - bioconda
      - defaults
      - conda-forge
    dependencies:
      - orthofinder=2.2.6

## Run pipeline in Johannesson's server

Get into the folder:

    $ cd /mnt/sda/johannesson_lab/podospora/4_PhylogenyPaper/1a_OrthoTrees

First, to get an idea of how the pipeline looks like we can make a rulegraph:
    
    <!-- $ conda install -c pkgs/main graphviz=2.40.1 -->
    $ snakemake --snakefile OrthoTrees.smk --configfile OrthoTrees_config.yaml --rulegraph | dot -Tpng > rulegraph.png

Run the pipeline:

    $ screen -R phylo
    $ conda activate LorePhylogenetics
    $ snakemake --snakefile OrthoTrees.smk --configfile OrthoTrees_config.yaml -p -j 35 --keep-going --use-conda &> OrthoTrees.log &


