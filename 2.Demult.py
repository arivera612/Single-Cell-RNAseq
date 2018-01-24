import sys,os
from sets import Set

# Function that calculates edit distance of 2 strings
def editDistance(seq1,seq2):
	distance = 0
	
	for i in range(0,min(len(seq1),len(seq2))):
		if seq1[i]!=seq2[i]:
			distance=distance+1
	return distance

# Function to stagger string into pieces
def stagger(seq,size):
	output = []
	for i in range(0,len(seq)-size+1):
		output.append(seq[i:i+size-1])
	return output

possibleBarcodes = ["AAAGAA","AACAGC","AACGTG","AAGCCA","AAGTAT","AATTGG","ACAAGG","ACCCAA","ACCTTC","ACGGAC","ACTGCA","AGACCC","AGATGT","AGCACG","AGGTTA","AGTAAA","AGTCTG","ATACTT","ATAGCG","ATATAC","ATCCGG","ATGAAG","ATTAGT","CAACCG","CAAGTC","CACCAC","CACTGT","CAGACT","CAGGAG","CATAGA","CCACGC","CCGATG","CCGTAA","CCTCTA","CGAAAG","CGAGCA","CGCATA","CGGCGT","CGGTCC","CGTTAT","CTAGGT","CTATTA","CTCAAT","CTGTGG","CTTACG","CTTGAA","GAAATA","GAAGGG","GACTCG","GAGCTT","GAGGCC","GAGTGA","GATCAA","GCCAGA","GCCGTT","GCGAAT","GCGCGG","GCTCCC","GCTGAG","GCTTGT","GGACGA","GGATTG","GGCCAT","GGGATC","GGTAGG","GGTGCT","GTACAG","GTCCTA","GTCGGC","GTGGTG","GTTAAC","GTTTCA","TAAGCT","TAATAG","TACCGA","TAGAGG","TATTTC","TCAGTG","TCATCA","TCCAAG","TCGCCT","TCGGGA","TCTAGC","TGAATT","TGAGAC","TGCGGT","TGCTAA","TGCTAA","TGGCAG","TGTGTA","TGTTCG","TTAAGA","TTCGCA","TTCGCA","TTCTTG","TTGCTC","TTGGAT","TTTGGG"]
#print possibleBarcodes[81]
# Input fastq files
print "Opening "+sys.argv[1]
read1File = open(sys.argv[1],'r')

print "Opening "+sys.argv[2]
read2File = open(sys.argv[2],'r')

sampleOutputPrefix = sys.argv[3]
outputprefix = sys.argv[4]
outputBarcodeFile = sys.argv[5]

sampleBarcodeFile = open(sys.argv[6],'r')

sampleBarcodes = []
Barcodes = {}
BarcodeUMIs = {}
Reads = {}
TotalBCs = {}
for line in sampleBarcodeFile:
	splitz = line.strip().split('\t')
	temp = [splitz[0],splitz[1]]
	sampleBarcodes.append(temp)
	Barcodes[splitz[1]]={}
	BarcodeUMIs[splitz[1]]={}
	Reads[splitz[1]]=[]
	TotalBCs[splitz[1]]=0
# For each read:
index = 0
read2=""
read1=""
read2storage=""
lineindex=0
read2Barcode = ""
skip = 0
for line in read1File:
	lineindex=lineindex+1
	if lineindex%100000==0:
		print lineindex
#	print index
#	print line
#	print read1
	line2 = read2File.readline()
	read2=read2+line2.strip()+"\n"
	if index == 1:
		read1 = line
	if index == 0:
		splitz = line.strip().split(':')
		read2Barcode = splitz[len(splitz)-1]
		found = 0
		for pair in sampleBarcodes:
			distance = editDistance(read2Barcode,pair[1])
			if distance <= 2:
				read2Barcode = pair[1]
				found = 1
				break
		if found == 0:
			skip = 1
	index=index+1
	if index == 4:
		read2storage=read2
		read2=""
		index = 0
	#	print read1
		if skip == 1:
			skip = 0
#			print "Read 2 Barcode not recognized: "+read2Barcode
			continue
		if len(read1) < 64:
			#print "Read too short"
			continue
# Stagger Read 1 into 15bp pieces
		staggerList = stagger(read1,15)
# Find Linker sequence 1 piece: TAGCCATCGCATTGC (Allowing 1 edit distance)
		Linker1Pos = -1
		for i in range(0,len(staggerList)):
			distance = editDistance(staggerList[i],"TAGCCATCGCATTGC")
			if distance <= 1:
				Linker1Pos = i
				break
# Make sure that Linker starts more than bp 6, continue if false
		if Linker1Pos < 6:
			#print "Linker 1 in wrong position: "+str(Linker1Pos)
			continue

# Find Linker sequence 2 piece: TACCTCTGAGCTGAA (Allowing 1 edit distance)
		Linker2Pos = -1
                for i in range(0,len(staggerList)):
                        distance = editDistance(staggerList[i],"TACCTCTGAGCTGAA")
                        if distance <= 1:
                                Linker2Pos = i
                                break

# Make sure that Linker starts are 21 bp apart, continue if false
		if Linker2Pos == -1 or Linker2Pos-21!=Linker1Pos:
			#print "Linker 1/2 in wrong position: "+str(Linker1Pos)+" "+str(Linker2Pos)
			continue
# Calculate BC positions
		BC1Pos = Linker1Pos-6
		BC2Pos = Linker1Pos+15
		BC3Pos = Linker1Pos+36
		ACGPos = Linker1Pos+42
		UMIPos = Linker1Pos+45
		GACTTPos = Linker1Pos+53
		EndPos = BC1Pos+64
		if len(read1) < EndPos:
			#print "Read ran off end"
			continue
# Make sure that ACG (allowing 1 edit distance) is after BC3, continue if false
		checkACGDistance = editDistance(read1[ACGPos:ACGPos+3],"ACG")
#		print "ACG "+read1[ACGPos:ACGPos+3]
		if checkACGDistance > 1:
			#print "ACG distance too far: "+str(checkACGDistance)+" "+read1[ACGPos:ACGPos+3]
			continue
# Make sure that GACTT (allowing 1 edit distance) is after UMI (8bp after ACG after BC3), continue if false
		checkGACTTDistance = editDistance(read1[GACTTPos:GACTTPos+5],"GACTT")
#		print "GACTT "+read1[GACTTPos:GACTTPos+5]
		if checkGACTTDistance > 1:
			#print "GACTT distance too far: "+str(checkGACTTDistance)+" "+read1[GACTTPos:GACTTPos+5]
			continue
# Calculate distance to each possible barcode block for BC1 (6 bp before Linker 1), BC2 (6bp after Linker1), BC3 (6bp after Linker2)
		BC1 = read1[BC1Pos:BC1Pos+6]
#		print "BC1 "+read1[BC1Pos:BC1Pos+6]
		barcode1Index = -1
		for i in range(0,len(possibleBarcodes)):
			checkBarcodeDistance = editDistance(BC1,possibleBarcodes[i])
			if checkBarcodeDistance <= 1:
				barcode1Index = i
				break
		if barcode1Index != -1:
			BC1 = possibleBarcodes[barcode1Index]
		else:
			continue
		BC2 = read1[BC2Pos:BC2Pos+6]
                barcode2Index = -1
                for i in range(0,len(possibleBarcodes)):
                        checkBarcodeDistance = editDistance(BC2,possibleBarcodes[i])
                        if checkBarcodeDistance <= 1:
                                barcode2Index = i
                                break
                if barcode2Index != -1:
                        BC2 = possibleBarcodes[barcode2Index]
		else:
			continue
		BC3 = read1[BC3Pos:BC3Pos+6]
		#print BC3
                barcode3Index = -1
                for i in range(0,len(possibleBarcodes)):
		#	print str(i)+" "+str(len(possibleBarcodes[i]))
                        checkBarcodeDistance = editDistance(BC3,possibleBarcodes[i])
                        if checkBarcodeDistance <= 1:
                                barcode3Index = i
                                break
                if barcode3Index != -1:
                        BC3 = possibleBarcodes[barcode3Index]
		else:
			continue

# If any barcode is more than 1 edit distance from it's correct form, leave it in

# Generate corrected barcode by combining corrected BC1, BC2, and BC3
		BC=BC1+BC2+BC3
		UMI = read1[UMIPos:UMIPos+8]

# Store UMIs for each barcode, if UMI repeated, continue
		if BC not in Barcodes[read2Barcode]:
			Barcodes[read2Barcode][BC]=TotalBCs[read2Barcode]
			TotalBCs[read2Barcode]=TotalBCs[read2Barcode]+1
			BarcodeUMIs[read2Barcode][BC]=Set()
			Reads[read2Barcode].append([])
		if UMI in BarcodeUMIs[read2Barcode][BC]:
			continue
		else:
			BarcodeUMIs[read2Barcode][BC].add(UMI)
# Output read 2 into barcode fastq file
		BarcodeIndex = Barcodes[read2Barcode][BC]
		Reads[read2Barcode][BarcodeIndex].append(read2storage)
for pair in sampleBarcodes:
	if not os.path.exists(sampleOutputPrefix+pair[0]):
		os.makedirs(sampleOutputPrefix+pair[0])
	for i in range(0,len(Reads[pair[1]])):
		outfile=open(sampleOutputPrefix+pair[0]+"/"+outputprefix+str(i)+".fastq",'w')
		for read in Reads[pair[1]][i]:
			outfile.write(read)
		outfile.close()

# Output cell ID -> Barcode table
for pair in sampleBarcodes:
	outfile = open(sampleOutputPrefix+pair[0]+"/"+outputBarcodeFile,'w')
	for BC in Barcodes[pair[1]]:
		outfile.write(outputprefix+str(Barcodes[pair[1]][BC])+'\t'+BC+'\n')
	outfile.close()
