# -*- snakemake -*-

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from ete3 import Tree

### TreeCertainty: Exploring the phylogenetic signal amongst orthologs of the *Podospora* complex
#############################################################################

# The goal of this pipeline is to assess the levels of conflict in the ML and
# coalescent topologies, using a number of metrics including the Internode
# Certainty (IC), Tree Certainty (TC) (Kobert et al. 2016 Mol. Biol. Evol.
# 33(6):1606–1617), a polytomy test, and the effect of using different degrees
# of phylogenetic signal in the individual genes.

#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/12/29
# ++++++++++++++++++++++++++++++++++++++++++++++
# Version 2

# -------------------------------------------------
AllSamples = config["SampleIDs"] # One per species! It already includes "Podan2" and "comataT"
samples = config["onerepperspp"] # One per species! It already includes "Podan2" and "comataT"
podan2gff = config["podan2gff"]

## --- Files from OrthoTrees.smk ---
# Output of OrthoFinder in csv
alignmentspath = config["alignmentspath"]
treespath = config["treespath"]
panames = config["panames"]
num_orthos = config["num_orthos"]

## --- Parameters ---
MINsupport = config["MINsupport"]

## --- Scripts ---
TreeCertainty = config["TreeCertainty"]
TC_env = config["TC_env"]
astral = config["astral"]

# -------------------------------------------------

def orthonums(num_orthos):
	allnums = []
	for n in range(1, num_orthos + 1):
		number = "{0:04d}".format(n)
		allnums.append(number)
	return allnums	


rule all:
	input:
		"results/TCdistributions.pdf",
		expand("results/TCgenes_top{MINrelTC}.txt", MINrelTC = [0.50, 0.70, 0.75]),
		"results/TCgenes_chrs.pdf",
		expand("concatenated/opinionatedgenes{MINrelTC}.treefile", MINrelTC = [0.50, 0.70, 0.75]), # concatenated Tree
		expand("astral/astral_top{MINrelTC}.bb{MINsupport}.poly.tre", MINrelTC = [0.50, 0.70, 0.75], MINsupport = MINsupport), # Astral Tree plus some extras
		expand("astral/astral_top{MINrelTC}.bb{MINsupport}.scored.tre", MINrelTC = [0.50, 0.70, 0.75], MINsupport = MINsupport), # Astral Tree plus some extras
		expand("results/RAxML_IC_Score_BranchLabels.opinionatedgenes{MINrelTC}.treefile", MINrelTC = [0.50, 0.70, 0.75], MINsupport = MINsupport), # Astral Tree plus some extras


rule reducealignments:
	""" Keep only one representative of each species, that is, keep only the ones in onerepperspp """
	input:
		alignmentspath + "/orthologs{n}.fas"
	output:
		"alignments/orthologs{n}_reduced.fas"
	run:
		# Read the fasta
		records_dict = SeqIO.to_dict(SeqIO.parse(input[0], "fasta", generic_dna))

		output_handle = open(output[0], "w")
		for strain in samples:
			# Get all sequences that match that strain
			matching = [seq for seq in records_dict.keys() if strain in seq]

			# Print the fasta sequences of those records
			for seq in matching:
				thiseq = records_dict[seq]
				thiseq.id = strain # Rename it
				thiseq.description = '' # Get read of annoying description
				SeqIO.write(thiseq, output_handle, "fasta")


rule IQTree: 
	""" Run IQTree for each ortholog group """
	input:
		"alignments/orthologs{n}_reduced.fas"
	output:
		"iqtree/orthologs{n}/orthologs{n}_reduced.treefile",
		"iqtree/orthologs{n}/orthologs{n}_reduced.ufboot",
	params:
		bootstraps = 1000, # UFBoot
		threads = 2
	shell:
		"""
		iqtree -s {input} -m MFP -seed 1234 -bb {params.bootstraps} -nt {params.threads} -bnni -pre "iqtree/orthologs{wildcards.n}/orthologs{wildcards.n}_reduced" -safe
		"""
		# -bnni to reduce the risk of overestimating branch supports with UFBoot due to severe model violations


rule InternodeCertainty_pergene:
	""" Calculate the internode Certainty of each gene tree """
	input:
		mltree = "iqtree/orthologs{n}/orthologs{n}_reduced.treefile",
		trees = "iqtree/orthologs{n}/orthologs{n}_reduced.ufboot",
	output:
		"IC/RAxML_info.orthologs{n}_reduced.treefile",
		"IC/RAxML_IC_Score_BranchLabels.orthologs{n}_reduced.treefile"
	shell:
		"raxmlHPC -f i -t {input.mltree} -z {input.trees} -m GTRCAT -w $PWD/IC -n orthologs{wildcards.n}_reduced.treefile"

def extractnum(line):
	g = re.findall(r'(.*): (.*)', line) # Notice the space after the colon
	return g[0][1]


rule extractTC:
	""" Make a table with the TC values """
	input:
		expand("IC/RAxML_info.orthologs{n}_reduced.treefile", n = orthonums(num_orthos))
	output:
		"results/TCgenes.txt"
	run: # I had to do it in python because the expansion of all the files kills the limit in bash
		header = "Ortholog\tTC\tRelTC\tTCA\tRelTCA\n"
		openoutput = open(output[0], 'w')
		openoutput.write(header)
	
		for file in input:
			openfile = open(file, 'r')

			ortho = re.findall(r'IC/RAxML_info.(.*)_reduced.treefile', file)[0]

			for line in openfile:
				if "Tree certainty for this tree:" in line:
					TC = extractnum(line)
				elif "Relative tree certainty for this tree:" in line:
					RelTC = extractnum(line)
				elif "Tree certainty including all conflicting bipartitions (TCA) for this tree:" in line:
					TCA = extractnum(line)
				elif "Relative tree certainty including all conflicting bipartitions (TCA) for this tree:" in line:
					RelTCA = extractnum(line)
				
			openoutput.write(f"{ortho}\t{TC}\t{RelTC}\t{TCA}\t{RelTCA}\n")

			openfile.close()


rule makegenetable:
	""" Make a table with information of the genes """
	input:
		podan2gff
	output:
		"results/GenesInfo.txt"
	shell:
		""" 
		grep -P '\\tgene' {input} | awk ' {{ attr=$9; id=gensub(/(.*);Name=([a-zA-Z0-9\_]*);(.*)/, "\\\\2", "g", attr); print $1,$4,$5,$7,id }} ' > {output}
		"""

rule plot_certainty:
	""" Plot the TC distribution and get a table of the best genes """
	input: # In this order
		"results/TCgenes.txt",
		"results/GenesInfo.txt",
		panames
	output: # In this order
		"results/TCdistributions.pdf",
		"results/TCgenes_top0.5.txt",
		"results/TCgenes_top0.7.txt",
		"results/TCgenes_top0.75.txt",
		"results/TCgenes_chrs.pdf",
	conda: 
		TC_env
	script:
		TreeCertainty

# From Salichos & Rokas 2013:
# "Internode-certainty values near zero indicate the presence of an almost
# equally supported bipartition that conflicts with the inferred internode,
# whereas values close to one indicate the absence of conflict."

# The difference between IC and ICA is that IC compares to the second most
# common bipartition, while ICA compares to all other bipartitions.


rule concatenation:
	""" Concatenate all ortholog alignments"""
	input:
		orthos = "results/TCgenes_top{MINrelTC}.txt",
	output:
		"concatenated/opinionatedgenes{MINrelTC}.fa"
	run: # similar to script fastaconcat.py
		# First find the alignments
		chosenorthos = []
		for line in open(input.orthos, 'r'):
			orthonumber = line.rstrip("\n")
			chosenorthos.append(orthonumber)
		# Make a list of the alignment file names
		alignments = [alignmentspath + f"/{i}.fas" for i in chosenorthos]

		# Read all fastas into lists with biopython objects
		allfastas = []
		for fastafile in alignments:
			thisfasta = [seq_record for seq_record in SeqIO.parse(fastafile, "fasta")]
			allfastas.append(thisfasta)

		# Make a dictionary of the sequences, starting with empty ones
		masternames = {} # A list of all names in all files
		for sample in AllSamples:
			masternames[sample] = Seq("")

		# Concatenate
		for lista in allfastas:
			thislen = len(lista[0].seq) # len of alignment in case I need to add missing data

			for name in masternames.keys():
				thisnameseq = "-"*thislen # assume it's not there

				for seq in lista:
					if name in seq.id: # it is there!
						thisnameseq = seq.seq
						break # you found it, so stop looping
				# Concat the sequence
				masternames[name] = masternames[name] + thisnameseq
		
		# Print concatenated sequence
		result = open(output[0], 'w')
		for sample in masternames.keys():
			result.write(">" + sample + '\n')
			result.write( str(masternames[sample]) + '\n')

rule IQTreeConcat:
	""" Make a tree of the concatenated alignment """
	input:
		"concatenated/opinionatedgenes{MINrelTC}.fa"
	output:
		"concatenated/opinionatedgenes{MINrelTC}.treefile"
	params:
		bootstraps = 1000, # UFBoot
		threads = 10,
	shell:
		"iqtree -s {input} -m MFP -seed 1234 -bb {params.bootstraps} -nt {params.threads} -bnni -pre 'concatenated/opinionatedgenes'{wildcards.MINrelTC}"

# ------ Trees for ASTRAL -------

rule getstrongtrees:
	""" Get the trees of the opinionated genes (from OrthoTrees.smk) """
	input:
		orthos = "results/TCgenes_top{MINrelTC}.txt",
	output:
		temp("astral/TCgenes_top{MINrelTC}.tre")
	params:
		treespath = treespath
	shell:
		"""
		for line in $(cat {input})
			do
				cat "{params.treespath}/$line/$line.fas.treefile" >> {output}
				# cat "{params.treespath}/$line/$line.treefile" >> {output}
			done
		"""

rule renametrees:
	""" Prepare trees for ASTRAL """
	input:
		"astral/TCgenes_top{MINrelTC}.tre"
	output:
		"astral/TCgenes_top{MINrelTC}.nom.tre"
	run:
		newtrees = []
		for line in open(input[0], 'r'):
			t = Tree(line)

			# Rename the leaves
			for node in t.iter_leaves():
				for sample in AllSamples:
					if sample in node.name:
						node.name = sample
						break
			# save
			newtrees += [t.write()]

		# # Print trees in a new file
		rootedtrees = open(output[0], 'w')
		for tree in newtrees:
			rootedtrees.write(tree + "\n")

# ---
rule InternodeCertainty_concat:
	""" Calculate the internode Certainty of the concatenated trees """
	input:
		concat = "concatenated/opinionatedgenes{MINrelTC}.treefile",
		trees = "astral/TCgenes_top{MINrelTC}.nom.tre"
	output:
		"results/RAxML_info.opinionatedgenes{MINrelTC}.treefile",
		"results/RAxML_IC_Score_BranchLabels.opinionatedgenes{MINrelTC}.treefile"
	shell:
		"raxmlHPC -f i -t {input.concat} -z {input.trees} -m GTRCAT -w $PWD/results -n opinionatedgenes{wildcards.MINrelTC}.treefile"
# ---

rule collapselowsupportbraches:
	""" Remove nodes in the tree that have low support using newick utilities """
	input:
		"astral/TCgenes_top{MINrelTC}.nom.tre"
	output:
		"astral/TCgenes_top{MINrelTC}.nom.bb{MINsupport}.tre"
	shell:
		"nw_ed {input} 'i & b<={wildcards.MINsupport}' o > {output}"
		# i matches internal nodes
		# b > 75 matches nodes whose label has a numerical value of 75 or more (if the label is numeric)
		# o   (splice Out) splice out node, and attach children to parent,
		# 		preserving branch lengths. This is useful for "opening" poorly
		#		supported nodes.

rule ASTRAL:
	""" Run ASTRAL on the collapsed trees """
	input:
		"astral/TCgenes_top{MINrelTC}.nom.bb{MINsupport}.tre"
	output:
		"astral/astral_top{MINrelTC}.bb{MINsupport}.tre"
	wildcard_constraints:
		MINsupport="\d+"
	params:
		astral = astral,
	shell:
		"java -jar {params.astral} -i {input} -o {output} 2> astral/astral_top{wildcards.MINrelTC}.bb{wildcards.MINsupport}.log"

rule scoring_astral:
	""" Produce extra branch support metrics for the Species tree of ASTRAL """
	input:
		genetrees = "astral/TCgenes_top{MINrelTC}.nom.bb{MINsupport}.tre",
		spptree = "astral/astral_top{MINrelTC}.bb{MINsupport}.tre"
	output:
		"astral/astral_top{MINrelTC}.bb{MINsupport}.scored.tre"
	wildcard_constraints:
		MINsupport="\d+"
	params:
		astral = astral,
	shell:
		"java -jar {params.astral} -q {input.spptree} -i {input.genetrees} -o {output} -t 8 2> astral/astral_top{wildcards.MINrelTC}.bb{wildcards.MINsupport}.scored.log"
		# -t 2 Full annotation (see https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md#extensive-branch-annotations)
		# -t 8 Alternative quartet topologies Outputs q1, q2, q3; these three values show quartet support (as defined in the description of -t 1) for the main topology, the first alternative, and the second alternative, respectively.
		#		Main topology: RL|SO, First alternative: RS|LO, and Second alternative: RO|LS

rule polytomy_test_astral:
	""" Test the null hypothesis of polytomi (see doi:10.3390/genes9030132) """
	input:
		genetrees = "astral/TCgenes_top{MINrelTC}.nom.bb{MINsupport}.tre",
		spptree = "astral/astral_top{MINrelTC}.bb{MINsupport}.tre"
	output:
		"astral/astral_top{MINrelTC}.bb{MINsupport}.poly.tre"
	wildcard_constraints:
		MINsupport="\d+"
	params:
		astral = astral,
	shell:
		"java -jar {params.astral} -q {input.spptree} -i {input.genetrees} -o {output} -t 10 2> astral/astral_top{wildcards.MINrelTC}.bb{wildcards.MINsupport}.poly.log"
		# -t 2 Full annotation (see https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md#extensive-branch-annotations)
