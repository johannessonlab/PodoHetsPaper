# -*- snakemake -*-

# import glob
from Bio import SeqIO
# from Bio import AlignIO
# from Bio.Alphabet import IUPAC, Gapped
# from Bio.Align import MultipleSeqAlignment
from ete3 import Tree

### PhyloHets: Making phylogenies of the het genes to detect ancestral polymorphism
#############################################################################

#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2020/01/26
# ++++++++++++++++++++++++++++++++++++++++++++++
# Version 1

# -------------------------------------------------
samples = config["SampleIDs"]
# Alignments
markers = config["markers"]
markerspath = config["markerspath"]

minbb = config["minbb"]

# Dotplot
pathsdotplotdata = config["pathsdotplotdata"]
seqsdotplots_ref = config["seqsdotplots_ref"]
seqsdotplots_query = config["seqsdotplots_query"]
gffs = config["gffs"]

# Scripts
dotplotenv = config["dotplotenv"]
dotplot = config["dotplot"]
# -------------------------------------------------

rule all:
	input:
		expand("results/{marker}.treefile", marker = markers),
		expand("mummer/{ref}-vs-{query}.coords", zip, ref = seqsdotplots_ref, query = seqsdotplots_query),
		"results/dotplot.pdf",
		"results/selfdotplot.pdf",


rule renamesamples:
	""" Rename sequences based on the strain """
	input:
		alignment = markerspath + "/{marker}"
	output:
		alignment = "alignments/{marker}"
	run:
		if 'hnwd' not in input.alignment:
			renamedseqs = [] # Setup an empty list

			for seq_record in SeqIO.parse(input.alignment, "fasta"):
				for sample in samples:
					if sample in seq_record.id:
						seq_record.id = sample
						seq_record.description = '' # otherwise it appends the old name, very annoying
						break # We found it, so let's move on
				renamedseqs.append(seq_record)

			# Write the sequences in the new fasta file
			SeqIO.write(renamedseqs, output.alignment, "fasta")
		else:
			shell("cp {input} {output}")


rule IQTree_marker:
	""" Make a tree of the marker """
	input:
		alignment = "alignments/{marker}"
	output:
		"iqtree/{marker}.bionj",
		"iqtree/{marker}.boottrees",
		"iqtree/{marker}.ckp.gz",
		"iqtree/{marker}.contree",
		"iqtree/{marker}.iqtree",
		# "iqtree/{marker}.log", # I probably need that in case it fails
		"iqtree/{marker}.mldist",
		"iqtree/{marker}.model.gz",
		"iqtree/{marker}.treefile"
	params:
		bootstraps = 100,
		threads = 4
	shell:
		"iqtree -s {input.alignment} -m MFP -seed 1234 -b {params.bootstraps} -nt {params.threads} -pre iqtree/{wildcards.marker} -redo"
		# -b <#replicates>     Bootstrap + ML tree + consensus tree (>=100)
		# -bb Ultrafast bootstraps (UFBoots)
		# -bnni to reduce the risk of overestimating branch supports with UFBoot due to severe model violations.


rule collect_results:
	input:
		"iqtree/{marker}.treefile"
	output:
		"results/{marker}.treefile"
	shell:
		"cp {input} {output}"


rule mummer:
	input:
		query = pathsdotplotdata + "{query}",
		reference = pathsdotplotdata + "{ref}",
	output:
		delta = "mummer/{ref}-vs-{query}.delta",
		deltafilter = "mummer/{ref}-vs-{query}.filter",
		coords = "mummer/{ref}-vs-{query}.coords",
		coordsfilter = "mummer/{ref}-vs-{query}.filter.coords",	
	params:
		threads = 3,
	shell:
		"""
		echo
		echo "MUMmer alignment ..."
		nucmer -b 200 -c 20 --maxmatch --nosimplify -p mummer/{wildcards.ref}-vs-{wildcards.query} {input.reference} {input.query} -t {params.threads}

		# Filter the delta
		delta-filter -q {output.delta} > {output.deltafilter}

		# To view a summary of all the alignments produced by NUCmer
		echo "Running show-coords"
		# For Ribbon http://genomeribbon.com/
		echo "...for Ribbon"
		show-coords -r -lTH {output.delta} > {output.coords}
		show-coords -r -lTH {output.deltafilter} > {output.coordsfilter}

		"""
		# --maxmatch      Use all anchor matches regardless of their uniqueness
		# --mum  Use anchor matches that are unique in both the reference and query
		# --mumreference  Use anchor matches that are unique in in the reference
        #           but not necessarily unique in the query (default behavior)
        # -c|mincluster   Sets the minimum length of a cluster of matches (default 65)
        # -b|breaklen     Set the distance an alignment extension will attempt to extend poor scoring regions before giving up (default 200)
        # --nosimplify    Don't simplify alignments by removing shadowed clusters. Use this option when aligning a sequence to itself to look for repeats

rule dotplot:
	""" Plot the TC distribution and get a table of the best genes """
	input: # In this order
		# Insertion
		expand("mummer/{ref}-vs-{query}.coords", ref = seqsdotplots_ref[0], query = seqsdotplots_query[0]),
		pathsdotplotdata + gffs[1],
		pathsdotplotdata + gffs[0],

		# Self-alignment 
		expand("mummer/{ref}-vs-{query}.coords", ref = seqsdotplots_ref[1], query = seqsdotplots_query[1]),
		pathsdotplotdata + gffs[2]

	output: # In this order
		"results/dotplot.pdf",
		"results/selfdotplot.pdf",
	conda: 
		dotplotenv
	script:
		dotplot
