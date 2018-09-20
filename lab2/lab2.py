"""
some useful website for Matplot lib:
https://matplotlib.org/users/pyplot_tutorial.html
https://matplotlib.org/examples/statistics/boxplot_demo.html
https://matplotlib.org/examples/statistics/histogram_demo_features.html

for BioPython:
https://biopython.org/wiki/SeqIO
https://biopython.org/wiki/Phylo
https://biopython.org/DIST/docs/api/Bio.Align.MultipleSeqAlignment-class.html

install matplotlib to visualize
"""

from Bio import SeqIO #input or output
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


SeqIO.parse(input_filename, format) #iterates over the contents of a file
SeqIO.write(sequences, output_handle, format) #writes a bunch of sequences to a file

#Seq is an object that contains a sequence—RNA, DNA, protein, etc.
#SeqRecord is an object that contains a Seq along with any metadata
#commonly supported format: fasta, fastq, genbank, stockholm, nexus, phylip



#example: FASTQ --> FASTA

my_recs = []
for seq_rec in SeqIO.parse("input.fastq", "fastq"):
	my_recs.append(seq_rec)
SeqIO.write(my_recs, open("output.fasta", "w"), "fasta")


#example: randomly generate a bunch of FASTA sequences

import random
my_recs = []
#generate 50 sequences
for seq_num in range(50):
	#each is 50 nucleotides long
	my_seq = ""
	for _ in range (50):
		my_seq += random.choice("ATCG")
	my_recs.append(SeqRecord(Seq = Seq(my_seq), id = "seq%d" % seq_num))
	SeqIO.write(my_recs, open("output.fasta", "w"), "fasta")

from Bio import AlignIO
alignment = AlignIO.read("seqs.aligned.fa", "fasta")
print ("Alignment length %i" % alignment.get_alignment_length())
for record in alignment:
	print(record.seq + " " + record.id) #return lists of alignments



import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
%matplotlib inline #to get """ipython""" to plot stuff in-line, this majic command needs to go somewhere

### Anatomy of a scatter plot ###
plot(x_axis_pts, y_axis_pts, format_string)
	#format_string is up to 3 characters consisting of (color) (point_style) (line_style)

axis([x_min, x_max, y_min, y_max])
	#show only a subset of the data by making range small
	#or put the data in context by making the range large

#example
%matplotlib inline
import matplotlib.pyplot as plt
plt.plot([1, 2, 3, 4], [1, 4, 9, 16], 'ro')
plt.axis([0, 6, 0, 20])
plot.show()


### Anatomy of a box plot ###
boxplot()

#example
%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
data = np.random.normal(size=(37,4))
line_out = plt.boxplot(data, labels=["A", "B", "C", "D"])

#for bar plot:
data = np.median(np.random.normal (loc=10, scale=20.0, size=(37,4)), axis=0)
	#generate random data and take the median of each set of 37 values
	#hieght on the graph equals to the median of each data set
line_out = plt.bar(height=data, x=["A", "B", "C", "D"])


### Quick tress in Newick form ###
Phylo.read(data, format)
#example
%matplotlib inline
from Bio import Phylo
from io import StringIO
tree = Phylo.read(StringIO ("A, (B, C), (D, E)"), "newick") #use StringIO passed in a file that handle out strings
Phylo.draw(tree)

### --------------------------------------------(end note)----------------------------------------------###
plt.xlabel('name')
plt.ylabel('name')
plt.title('title')
plt.bar()

### ----------------

from Bio import Phylo
tree = Phylo.read("filename.nwk", "newick") #file name and format
tree.ladderize()
Phylo.draw(tree) #tree shown in png
Phylo.draw_ascii(tree) #tree shown in terminal

from Bio import SeqIO
input_iter = SeqIO.parse("seqs.fa", "fasta")
#print("hu.32".seq)
for record in input_iter:
	if record.id == "hu.14" or record.id == "hu.31" or record.id == "hu.32": #limit the options
		print(record.seq)

from Bio import AlignIO
alignment = AlignIO.read("seqs.aligned.fa", "fasta")
print("Alignment length %i" % alignment.get_alignment_length()) #alignment length
for record in alignment:
     print(record.seq + " " + record.id)


import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

input_iter = SeqIO.parse("seqs.fa", "fasta")
for record in input_iter:
	if record.id == "hu.14" or record.id == "hu.31" or record.id == "hu.32": #limit the options
		print(record.seq.count("a"))

#---------bar plot----------#

from Bio import SeqIO
import numpy as np
import matplotlib as mlab
import matplotlib.pyplot as plt
input_iter = SeqIO.parse("seqs.aligned.fa", "fasta")
chosen_specimen = ["hu.39", "rh.58", "pi.2", "hu.17", "rh.40", "hu.14", 
                   "hu.43", "cy.2", "rh.35", "cy.3" ]
result = []
all_seq = []
for record in input_iter:
	if record.id in chosen_specimen:
		temp = record.seq
		all_seq.append(temp) #add sequence (as string) to a created list.
seq_len = len(all_seq[0])
for i in range(seq_len):
	temp = ""
	for seq in all_seq:
		temp += seq[i] #turn all the i-th nucleotide of all sequences to string
	result.append(temp)

percentAT = [(_.count("A") + _.count("T"))/len(chosen_specimen)for _ in result] # %AT_content of all nt
percentGC = [1 - _ for _ in percentAT]
plt.bar(np.arange(seq_len), percentAT)
plt.bar(np.arange(seq_len), percentGC, bottom=percentAT)
plt.xlabel("Position in Sequence")
plt.ylabel("% AT or % GC")
plt.xticks(np.arange(seq_len)) #add ,(text) if set text labels
plt.yticks(np.arange(1, step=0.1))
plt.show()

#----------box plot ----------#
cluster1 = []
cluster2 = []
cluster3 = []
cluster4 = []
input_iter = SeqIO.parse("seqs.fa", "fasta")
for record in input_iter:
	if record.id in ["hu.39", "rh.50", "rh.49","rh.57", "rh.51", 
					"rh.53", "rh.64", "rh.52", "rh.61", "rh.58"]:
		cluster1.append(len(record.seq))
	elif record.id in ["pi.2", "pi.3", "pi.1", "hu.17", "hu.6", 
						"bb.1", "bb.2", "rh.2", "rh.40", "hu.67",
						"hu.37", "hu.40", "hu.66", "hu.42", "hu.41", "rh.38"]:
		cluster2.append(len(record.seq))
	elif record.id in ["rh.43", "hu.14", "hu.31", "hu.32", "hu.43", "hu.48",
						"hu.44", "hu.46", "cy.2", "rh.54", "rh.55", "rh.48", "rh.62"]:
		cluster3.append(len(record.seq))
	elif record.id in ["cy.2", "rh.54", "rh.55", "rh.48", "rh.62", "rh.35", "rh.36", "rh.37",
						"cy.3", "cy.6", "cy.4", "cy.5", "rh.13"]:
		cluster4.append(len(record.seq))
data = [cluster2, cluster2, cluster3, cluster4]
plt.boxplot(data, labels=["clust 1", "clust 2", "clust 3", "clust 4"], autorange=True, meanline=True)
plt.ylabel('Sequence length')
plt.show()




