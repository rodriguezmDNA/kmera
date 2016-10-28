#!/usr/bin/python
import re # Module for regular expressions
import time

# -- Start timer, open log file
start_time = time.time()

logFile = open('Motif.log', 'w')
print >> logFile, 'FASTA reading and kmer finding \n'  # or f.write('...\n')
#logFile.write('This is my log file \n')
#logFile.write('...\n')


# -- Where functions live
def slidingWindow(header,sequence, kmerSize, SeqLen,kmerDictionary):
	SeqStart = 0
	SeqEnd = (SeqLen - kmerSize) - 1 #-1 because of Python's 0-Index
	# Sliding window through sequence,
	if kmerSize <= SeqLen:
		print >> logFile, header, "ok"  # Write to log
		while SeqStart <= SeqEnd:
			kmer = sequence[SeqStart:SeqStart+kmerSize]
			#print (kmer)
			# -- Call function
			kmerCheck(kmerDictionary,kmer)
			# -*
			# Slide to next base:
			SeqStart += 1
	else:
		print >> logFile, header, "Kmer longer than sequence"  # Write to log
		
	return(kmerDictionary)


def kmerCheck(kmerDictionary,kmer):
	# Check kmer in dictionary. If found, increase count by 1. If not, create with starting count 1.
	if kmer in kmerDictionary.keys():
		#print ("True")
		kmerDictionary[kmer] += 1
	else:
		kmerDictionary[kmer] = 1
	return (kmerDictionary)



def readFASTA(fastafile,kmerSize,kmerDictionary):
	# Read file
	FASTA = open(fastafile)
	
	SeqCount = 0
	TotalLengths = 0
	for eachLine in FASTA:
		# Check FASTA header
		if eachLine[0] == '>':
				header = re.sub('\:$', '', str.split(eachLine)[1]) #Split the header, get only the AGI
				print (header)
				# Grab sequence
				sequence = FASTA.next() #Save sequence
				SeqLen = len(sequence)
				#
				SeqCount += 1 # Count each sequence read.
		
		TotalLengths +=  SeqLen #Sum lengths of sequences
		#print (sequence)
		# -- Call sliding window function
		slidingWindow(header,sequence, kmerSize, SeqLen,kmerDictionary)
		# -*
	
	#
	avgLength = TotalLengths/SeqCount
	print >> logFile, "\n--"
	print >> logFile, "Average sequence length:", avgLength
	print >> logFile, "K-mer size:", kmerSize
	# --
	return(kmerDictionary)
	FASTA.close()
	

def writeResults(kmerDictionary,resultFileName):
	# --
	# Write results to file as table
	# Open file
	resultsFile = open(resultFileName, 'w')
	resultsFile.write('Kmer \t Counts\n')

	for key in kmerDictionary.keys():
		resultsFile.write (key + "\t" + str(kmerDictionary[key]) + "\n")
	resultsFile.close()
	# -*


## -- Set initial parameters and variables
# Parameters
kmerSize = 8
fastafile = "../promoters_onlyNLP_NUE_All_v2.fa"
# Variables
kmerDictionary = {}
##

# Make names
ext = ".txt"
fasta_basename = re.sub('.fa|../', '', fastafile)
resultFileName = "kmers_" + str(kmerSize) + "_" + fasta_basename + ext


# Read FASTA
print >> logFile, "Processing:",fastafile  # Write to log
readFASTA (fastafile,kmerSize,kmerDictionary)
#print (kmerDictionary)

# Write results to file
print >> logFile, "Writing results to:",resultFileName  # Write to log
writeResults(kmerDictionary,resultFileName)



# -- Close files
logFile.close()
print ("\n")
print("--- Time taken: %s seconds ---" % (time.time() - start_time))

