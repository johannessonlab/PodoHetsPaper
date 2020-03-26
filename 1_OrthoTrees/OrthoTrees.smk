# -*- snakemake -*-

# from pytools.persistent_dict import PersistentDict
from glob import glob
from Bio import SeqIO
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment
from ete3 import Tree

### OrthoTrees: Getting orthologs of the Podospora complex
#############################################################################

# A general pipeline to get orthologous genes from the *Podospora* complex.
# Due to the unpredictable number of ortholog groups, the pipeline relies on a
# [checkpoint function]. I ran into a bug and I had to do some very awkward
# work around. 

# The idea is that I first get the orthogroups with OrthoFinder, and then, to
# avoid annotation issues, I BLAST back the Podan2 genes to all the other
# species. Thus, I retrieve the nucleotide sequences with introns. I can do
# this because the species are very closely related.

#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/12/19 - 25
# ++++++++++++++++++++++++++++++++++++++++++++++
# Version 3

# -------------------------------------------------
samples = config["SampleIDs"]
assembliespath = config["assembliespath"]

# References
podan2 = config["podan2"]
podan2gff = config["podan2gff"]
podan2genes = config["podan2genes"]
comataT = config["comataT"]
comataTgff = config["comataTgff"]

# Outgroup
outgroup = config["outgroup"]


# Scripts
gff2fasta = config["gff2fasta"]
orthogrs_parser = config["orthogrs_parser"]
query2hitseq = config["query2hitseq"]
fastaconcat = config["fastaconcat"]

AllSamples = samples + ["Podan2", "comataT"]
# -------------------------------------------------




def roottree(tree, outgroup, f = 0):
	t = Tree(tree, format = f)
	if len(outgroup) == 1:
		t.set_outgroup(outgroup[0])
	else:
		ancestor = t.get_common_ancestor(outgroup)
		t.set_outgroup(ancestor)
	return(t)

rule all:
	input:
		"results/SingleGeneTrees.tre"


# ------- PREPARE ALL DATA --------
## Make symlink of the genomes to work more easily
rule getgenomes:
	""" Make links to the assemblies """
	input: 
		assembliespath + "/{sample}.nice.fa",
	output:
		"genomes/{sample}.fa"
	shell:
		"ln -sf {input} {output}"

rule getrefgenomes:
	""" Make links to the assemblies """
	input: 
		podan2 = podan2,
		comataT = comataT,
	output:
		podan2 = "genomes/Podan2.fa",
		comataT = "genomes/comataT.fa",
	shell:
		"cat {input.podan2} | sed 's;>;>Podan2_;' > {output.podan2};"
		"cat {input.comataT} | sed 's;>;>comataT_;' > {output.comataT};"

# ----------------------------------
rule getprotsamples:
	""" Get CDS sequences """ 
	input:
		genome = assembliespath + "/{sample}.nice.fa",
		gff = lambda wildcards: glob(assembliespath + "/{sample}.nice*.gff3".format(sample = wildcards.sample)) # Dirty trick to expand * like in bash
	output:
		prots = "proteins/{sample}.fas",
	params:
		gff2fasta = gff2fasta
	shell:
		"""
		# Get CDS translated
		python {params.gff2fasta} {input.genome} {input.gff} --output proteins/{wildcards.sample} -t CDS -p -j --onlynames
	
		# Put the sample ID in the name of the sequences
		sed -i 's;>\\(.*\\);>\\1_{wildcards.sample};' {output.prots} 
		"""

rule getprotpodan:
	""" Get CDS sequences of reference genomes""" 
	input:
		genome = podan2,
		gff = podan2gff,
	output:
		prots = "proteins/Podan2.fas",
	params:
		gff2fasta = gff2fasta
	shell:
		"""
		# Get CDS translated
		python {params.gff2fasta} {input.genome} {input.gff} --output proteins/Podan2 -t CDS -p -j --onlynames
	
		# Put the sample ID in the name of the sequences
		sed -i 's;>\\(.*\\);>\\1_Podan2;' {output.prots} 
		"""

rule getprotcomata:
	""" Get CDS sequences of reference genomes""" 
	input:
		genome = comataT,
		gff = comataTgff,
	output:
		prots = "proteins/comataT.fas",
	params:
		gff2fasta = gff2fasta
	shell:
		"""
		# Get CDS translated
		python {params.gff2fasta} {input.genome} {input.gff} --output proteins/comataT -t CDS -p -j --onlynames
	
		# Put the sample ID in the name of the sequences
		sed -i 's;>\\(.*\\);>\\1_comataT;' {output.prots} 
		"""

# -----------------------------------------


rule orthofinder:
	""" Run OrthoFinder """
	# https://github.com/davidemms/OrthoFinder/blob/master/OrthoFinder-manual.pdf
	input:
		expand("proteins/{sample}.fas", sample = AllSamples)
	output:
		"orthofinder/Orthogroups.csv"
	params:
		threads = 30,
	conda: 
		"envs/orthofinder.yaml"
	shell:
		"""
		# Run Orthofinder
		orthofinder -f proteins -t {params.threads}

		# Move the new folder and change name
		mv proteins/Results_*/* orthofinder/
		"""

rule parseorthogroups1n:
	""" Parse the output of orthogroups """
	input:
		"orthofinder/Orthogroups.csv"
	output:
		podan1n = "filtering/Podan2_1n.txt",
		podan1nClean = "filtering/Podan2_1n.clean.txt",
	params:
		orthogrs_parser = orthogrs_parser,
		threads = 1,
	shell:
		"""
		# Filter the Orthogroups.csv for groups of one-to-one orthologs present in all samples
		{params.orthogrs_parser} {input} -n1 -b -o filtering
		
		# Clean the output a bit
		sed -i 's;_Podan2;;g' {output.podan1n}

		# Remove genes that are not starting with "Pa_"
		grep '^Pa' {output.podan1n} > {output.podan1nClean}
		"""


checkpoint orthofolders:
	""" Make dummy files for every orthogroup """
	input:
		orthologs = "filtering/Podan2_1n.clean.txt",
	output:
		directory("fastahits")
	run:
		shell("mkdir -p fastahits")
		shell("mkdir -p fastahits/dummies")
		num_lines = sum(1 for line in open(input[0])) # How many lines in the file?

		for n in range(1, num_lines + 1):
			number = "{0:04d}".format(n)
			fname = "fastahits/dummies/orthologs%s.dummy" % (number)
			tabopen = open(fname, 'w')

def dummiesoutput(wildcards):
	checkpoint_output = checkpoints.orthofolders.get(**wildcards).output[0] # "fastahits"

	return expand("fastahits/dummies/orthologs{i}.dummy",
		i=glob_wildcards(os.path.join(checkpoint_output, 'dummies/orthologs{i}.dummy')).i)

def iqtreeoutput(wildcards):
	checkpoint_output = checkpoints.orthofolders.get(**wildcards)#.output[0] # I actually don't need anything from this, just to confirm it's a checkpoint

	# Get the names of all the final alignments
	file = glob("filtering/Podan2_1n.clean.txt")[0]
	num_lines = sum(1 for line in open(file)) # How many lines in the file?
	finalnames = []
	for n in range(1, num_lines + 1):
		number = "{0:04d}".format(n)
		fname = "iqtree/orthologs%s/orthologs%s.fas.treefile" % (number, number)
		finalnames.append(fname)
	return finalnames

rule query2hitseq_persample:
	""" Get fasta files of each ortholog and for each sample"""
	input:
		dummy = expand("fastahits/dummies/orthologs{i}.dummy", i = lambda wildcards: glob_wildcards(os.path.join(checkpoints.orthofolders.get(**wildcards).output[0], 'dummies/orthologs{i}.dummy')).i ),-
		genome = "genomes/{sample}.fa",
		orthologs = "filtering/Podan2_1n.clean.txt",
		# dummy = expand("fastahits/dummies/orthologs{i}.dummy", i = ["{0:04d}".format(i) for i in range(1, sum([1 for line in open("filtering/Podan2_1n.clean.txt", 'r') ]) + 1 )] ),
	output:
		# expand("fastahits/orthologs{i}/{{sample}}.fa", i = lambda wildcards: glob_wildcards(os.path.join(checkpoints.orthofolders.get(**wildcards).output[0], 'dummies/orthologs{i}.dummy')).i )
		expand("fastahits/orthologs{i}/{{sample}}.fa", i = ["{0:04d}".format(i) for i in range(1, sum([1 for line in open("filtering/Podan2_1n.clean.txt", 'r') ]) + 1 )] )
	params: 
		query2hitseq = query2hitseq,
		refgenes = podan2genes,
	run:
		# Read the ortholog groups file
		tabopen = open(input.orthologs, 'r')
		tabs = [line.rstrip("\n") for line in tabopen] 			# Read tab file into a list
		
		# Get fasta
		count = 1
		for ortho in tabs:
			number = "{0:04d}".format(count)
			cmd1 = f"python {{params.query2hitseq}} {{params.refgenes}} {{input.genome}} -i {{ortho}} --tophit --extrabp 0 --temp fastahits/orthologs{{number}}/{{wildcards.sample}}/ >> fastahits/orthologs{{number}}/{{wildcards.sample}}.fa"
			shell(cmd1)

			count += 1

rule cathits:
	""" Put the sequences of all samples into a single fasta """
	input:
		# cathitsinput
		expand("fastahits/orthologs{{i}}/{sample}.fa", sample = AllSamples) # Notice only the sample wildcards gets expanded
	output:
		"fastahits/orthologs{i}.fas"
	shell:
		"cat {input} > {output}"


rule mafft:
	""" Align ortholog groups with MAFFT """
	input:
		ortholog = "fastahits/orthologs{i}.fas",
	output:
		"alignments/orthologs{i}.fas",
	params:
		threads = 6,
	shell:
		""" 
		mafft --thread {params.threads} --threadit 0 --adjustdirection --anysymbol --maxiterate 1000 --retree 1 --localpair {input.ortholog} > {output}
		"""


rule IQTreePerGene: # Careful, some trees don't have all the samples
	""" Run IQTree for each ortholog group """
	input:
		"alignments/orthologs{i}.fas"
	output:
		"iqtree/orthologs{i}/orthologs{i}.fas.treefile"
	params:
		bootstraps = 100,
		threads = 3
	shell:
		"""
		cd iqtree/orthologs{wildcards.i}
		ln -sf ../../{input} orthologs{wildcards.i}.fas
		iqtree -s orthologs{wildcards.i}.fas -m MFP -seed 1234 -b {params.bootstraps} -nt {params.threads} -bnni
		"""
		# -bnni to reduce the risk of overestimating branch supports with UFBoot due to severe model violations
		# -bb UFboot
		# -b standard nonparametric bootstrap

rule aggregate_trees:
	""" Collect all trees in a single file """ 
	# With this rule I trigger the formation of all other checkpoint files
	input:
		iqtreeoutput,
	output:
		"results/SingleGeneTrees.tre"
	# shell:
	# 	"cat {input} > {output}"
	run: # I had to do it like this because the line becomes so long that i get an error ([Errno 7] Argument list too long)
		# Save all trees in a single file
		alltrees = []
		for tree in input:
			t = Tree(tree)
			alltrees += [t.write()] # Save the newick line as is

		# Save all trees in a single file
		result = open(output[0], 'w')
		for tree in alltrees:
			result.write(tree + "\n")

