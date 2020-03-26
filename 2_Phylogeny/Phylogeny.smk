# -*- snakemake -*-

from ete3 import Tree
from Bio import SeqIO
from Bio.Seq import Seq

### Phylogeny: Making a phylogeny of the Podospora complex
#############################################################################

# A pipeline to make a phylogeny of *Podospora* complex. 

#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/12/25
# ++++++++++++++++++++++++++++++++++++++++++++++
# Version 1

# -------------------------------------------------
samples = config["SampleIDs"]

## --- Files from OrthoTrees.smk ---
alignmentspath = config["alignmentspath"]
alltrees = config["alltrees"]

# Output of OrthoFinder in csv
orthogroups = config["orthogroups"]

# Orthogroup equivalence in Podan codes
pa_orthogroups = config["pa_orthogroups"]
## ---------
# Outgroup
outgroup = config["outgroup"]
# Reference genes
podan2genes = config["podan2genes"]

# ----------
# Number of sample orthologs
SAMPLEsize = config["SAMPLEsize"]
# Support filter for collapsing branches
MINsupport= config["MINsupport"]

# ----------
# Scripts
orthogrs_parser = config["orthogrs_parser"]
query2hitseq = config["query2hitseq"]

astral = config["astral"]
mappingfile = config["mappingfile"]
# ----------

AllSamples = samples
# -------------------------------------------------

rule all:
	input:
		"results/RAxML_info.allonetoones.tre",
		"results/alltrees.nom.0intraspp.tre",
		expand("astral/astral.bb{MINsupport}.scored.tre", MINsupport = MINsupport), # For ASTRAL
		expand("astral/astral.bb{MINsupport}.poly.tre", MINsupport = MINsupport), # For ASTRAL
		"concatenated/allonetoones.tre", # Concatenation tree

# ------ Concatenated tree ------

rule parseorthogroups1n:
	""" Parse the output of orthogroups """
	input:
		orthogroups
	output:
		expand("filtering/Podan2_1n.sample{SAMPLEsize}.txt", SAMPLEsize = SAMPLEsize)
	params:
		orthogrs_parser = orthogrs_parser,
		threads = 1,
		SAMPLEsize = SAMPLEsize
	shell:
		"""
		# Filter the Orthogroups.csv for groups of one-to-one orthologs present in all samples
		{params.orthogrs_parser} {input} -n1 -b -o filtering -s {params.SAMPLEsize}
		
		# Clean the output a bit
		sed -i 's;_Podan2;;g' filtering/Podan2_1n.txt

		# Change name
		mv filtering/Podan2_1n.txt {output}
		"""

# Make a dictionary of all Pa genes and the alignment number
orthodict = {}
orthofile = open(pa_orthogroups, 'r')
n = 1
for line in orthofile:
	panumber = line.rstrip("\n")
	orthonumber = "{0:04d}".format(n)
	orthodict[panumber] = orthonumber
	n += 1
orthofile.close()

def roottree(tree, outgroup, f = 0):
	t = Tree(tree, format = f)
	if len(outgroup) == 1:
		t.set_outgroup(outgroup[0])
	else:
		ancestor = t.get_common_ancestor(outgroup)
		t.set_outgroup(ancestor)
	return(t)

rule concatenation:
	""" Concatenate all the one-to-one ortholog alignments"""
	input:
		orthos = expand("filtering/Podan2_1n.sample{SAMPLEsize}.txt", SAMPLEsize = SAMPLEsize),
	output:
		"concatenated/allonetoones.fa"
	run: # similar to script fastaconcat.py
		# First find the alignments
		chosenumbers = []
		for line in open(input.orthos[0], 'r'):
			panumber = line.rstrip("\n")
			chosenumbers.append(orthodict[panumber])
		# Make a list of the tree file names
		alignments = [alignmentspath + f"/orthologs{i}.fas" for i in chosenumbers]

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
		"concatenated/allonetoones.fa"
	output:
		"concatenated/allonetoones.tre"
	params:
		bootstraps = 1000, # UFBoot
		threads = 10
	shell:
		"""
		iqtree -s {input} -m MFP -seed 1234 -bb {params.bootstraps} -nt {params.threads} -bnni -pre "concatenated/allonetoones"
		mv concatenated/allonetoones.treefile {output}
		"""

# ------ Trees for ASTRAL -------

rule renametrees:
	""" Prepare trees for ASTRAL """
	input:
		alltrees
	output:
		"results/alltrees.nom.tre"
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

		# Print trees in a new file
		rootedtrees = open(output[0], 'w')
		for tree in newtrees:
			rootedtrees.write(tree + "\n")

rule collapselowsupportbraches:
	""" Remove nodes in the tree that have low support using newick utilities """
	input: 
		"results/alltrees.nom.tre"
	output:
		"astral/alltrees.nom.bb{MINsupport}.tre"
	params:
		minsupport = MINsupport
	shell:
		"nw_ed {input} 'i & b<={params.minsupport}' o > {output}"
		# i matches internal nodes
		# b > 75 matches nodes whose label has a numerical value of 75 or more (if the label is numeric)
		# o   (splice Out) splice out node, and attach children to parent,
		# 		preserving branch lengths. This is useful for "opening" poorly
		#		supported nodes.

rule ASTRAL:
	""" Run ASTRAL on the collapsed trees """
	input:
		"astral/alltrees.nom.bb{MINsupport}.tre"
	output:
		"astral/astral.bb{MINsupport}.tre"
	wildcard_constraints:
		MINsupport="\d+"
	params:
		astral = astral,
		# mappingfile = mappingfile, # what sample belongs to what species
		log = "astral/astral.bb" + str(MINsupport) + ".log",
	shell:
		"java -jar {params.astral} -i {input} -o {output} 2> {params.log}"
		# "java -jar {params.astral} -i {input} -a {params.mappingfile} -o {output} 2> {params.log}"

rule scoring_astral:
	""" Produce extra branch support metrics for the Species tree of ASTRAL """
	input:
		genetrees = "astral/alltrees.nom.bb{MINsupport}.tre",
		spptree = "astral/astral.bb{MINsupport}.tre"
	output:
		"astral/astral.bb{MINsupport}.scored.tre"
	wildcard_constraints:
		MINsupport="\d+"
	params:
		astral = astral,
		log = "astral/astral.bb" + str(MINsupport) + ".scored.log",
	shell:
		"java -jar {params.astral} -q {input.spptree} -i {input.genetrees} -o {output} -t 8 2> {params.log}"
		# -t 2 Full annotation (see https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md#extensive-branch-annotations)
		# -t 8 Alternative quartet topologies Outputs q1, q2, q3; these three values show quartet support (as defined in the description of -t 1) for the main topology, the first alternative, and the second alternative, respectively.
		#		Main topology: RL|SO, First alternative: RS|LO, and Second alternative: RO|LS

rule polytomy_test_astral:
	""" Test the null hypothesis of polytomi (see doi:10.3390/genes9030132) """
	input:
		genetrees = "astral/alltrees.nom.bb{MINsupport}.tre",
		spptree = "astral/astral.bb{MINsupport}.tre"
	output:
		"astral/astral.bb{MINsupport}.poly.tre"
	wildcard_constraints:
		MINsupport="\d+"
	params:
		astral = astral,
		log = "astral/astral.bb" + str(MINsupport) + ".poly.log",
	shell:
		"java -jar {params.astral} -q {input.spptree} -i {input.genetrees} -o {output} -t 10 2> {params.log}"
		# -t 2 Full annotation (see https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md#extensive-branch-annotations)

# ------ Opinionated trees ------

## Species definitions, this could be read from the spp_mapping.txt file but I'm being lazy.

def search_by_size(node, size):
    "Finds nodes with a given number of leaves"
    matches = []
    for n in node.traverse():
       if len(n) == size:
          matches.append(n)
    return matches

rule cleanIntraSppSupport:
	""" Replace the intra-spp support branches with 0 """
	input:
		map = mappingfile,
		trees = "results/alltrees.nom.tre",
	output:
		"results/alltrees.nom.0intraspp.tre"
	run:
		# Turn the ASTRAL mapping file into a dictionary
		mappingfile = open(input.map, 'r')
		multiindivspp = {}
		tabs = [line.rstrip("\n").split(":") for line in mappingfile] # Read every line in mapping file

		for tab in tabs:
			sp = tab[0]
			strains = [s for s in tab[1].split(",")]
			if len(strains) > 1: multiindivspp[sp] = strains # Keep only species with more than one individual

		# Modify the trees
		newtrees = []
		for line in open(input.trees, 'r'):
			t = Tree(line)
			t.set_outgroup(outgroup[0]) # Root the tree so I get anserina monophyletic
			# print(search_by_size(t, size=6))

			# for node in t.iter_descendants("postorder"): # Iterate over a tree excluding the root node
			for node in t.traverse(): # Iterate over a tree 
				if not node.is_leaf():
					descendants = [c.name for c in node.iter_leaves()] # What are the strains that come from this node?
					for sp in multiindivspp.keys():
						if set(descendants).issubset(set(multiindivspp[sp])): # Are the descendants a subset of the species?
							node.support = 0

			# save
			newtrees += [t.write(format = 0)]

		# Print trees in a new file
		rootedtrees = open(output[0], 'w')
		for tree in newtrees:
			rootedtrees.write(tree + "\n")	


rule InternodeCertainty:
	""" Calculate the internode Certainty of the ML and ASTRAL trees """
	input:
		concat = "concatenated/allonetoones.tre",
		trees = "results/alltrees.nom.tre"
	output:
		"results/RAxML_info.allonetoones.tre",
		"results/RAxML_IC_Score_BranchLabels.allonetoones.tre"
	shell:
		"raxmlHPC -f i -t {input.concat} -z {input.trees} -m GTRCAT -w $PWD/results -n allonetoones.tre"

