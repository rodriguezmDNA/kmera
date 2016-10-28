#!/usr/bin/python
import re # Module for regular expressions
import time

f = open('Motif.log', 'w')
print >> f, 'This is my log file \n'  # or f.write('...\n')
#f.write('This is my log file \n')
#f.write('...\n')
f.close()


start_time = time.time()



# Read file
promoters_file = open("../promoters_onlyNLP_NUE_All_v2.fa")

## Some info about regexes in Python
# <re.match> checks for a match only at the beginning of the string, 
# <re.search> checks anywhere in the string (Perl does by default).
# <re.findall> finds all regex matches


# Read the fasta file
def fn_Motif(FastaFile,regexpr):
	totalMatches = 0
	avgSeqLen = 0
	## Where the magic starts
	for line in FastaFile:
		#print header[0]
		sequence = FastaFile.next()
		# Check FASTA header
		if line[0] == '>':
			header = re.sub('\:$', '', str.split(line)[1]) #Split the header, get only the AGI
			print (header) 
		# Do stuff here:	
		#targetGenes.append(header)
		MotifFinder = re.findall(regexpr,sequence,re.I)
		match_times = len(MotifFinder)
		seqLen = len(sequence)
		print (MotifFinder)
		allMotifs[header] = (MotifFinder)
		# Add stuff to the dictionary
		targetGenes[header]=(match_times,seqLen)
		totalMatches += match_times
		avgSeqLen += seqLen
		#print (line)
		return (targetGenes)
		return (totalMatches)
		return (avgSeqLen)

targetGenes = {} # Use dictionaries to be able to use strings as indexes
motifs_found = ()
allMotifs = {}

targetGenes.keys()

motif="a.{9,10}t{3,7}" #Try with match
fn_Motif(promoters_file,motif)

#avgSeqLen = avgSeqLen/len(targetGenes)
#print (avgSeqLen)


"""
avgSeqLen = avgSeqLen/len(targetGenes)


#motif = ".GAC.CTT.{0,10}AAG" # NCTTNNNNNNNNNNAAG
#motif="AATT" #Try with match


fn_Motif(promoters_file,motif)

#print (MotifFinder)
print ("\n")
print ("Average gene length:" + str(avgSeqLen))
print ("Total number of matches:" + str(totalMatches))

#outFile = open('targetGenes.txt', 'w')

with open ("targetGenes.tsv", "w") as f:
    f.write("" + "\t" + "match_times" + "\t" + "SeqLength" + "\n")
    for key, value in targetGenes.items():
        f.write("{}\t{}\n".format(key, "\t".join([str(i) for i in value])))


#




print ("\n")
print("--- %s seconds ---" % (time.time() - start_time))


## -- match function
#m = re.match(motif,sequence,re.I) # Search for the motif
## -- *

## -- findall creates a list
#m = re.findall(motif,sequence,re.I)
#print (m)

## -- *

# -- Iteratively
#m = re.finditer(motif,sequence,re.I)
#for each in m:
#	print each.group(0)
# -- *

#
#print (m.groups)




# -- 
#m = re.search('(?<=abc)def', 'abcdef')
# --


"""
