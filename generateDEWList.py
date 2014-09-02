#!/usr/bin/env python


'''
DEW: Differentially Expressed Words
    FUNCTION HELP:	N: get DEW List
    IN:			Class A (multi)fasta files
			Class B (multi)fasta files
    OUT: 		ods doc with sorted DEW

    TIME COMPLEXITY: 	O(boy)
			Probably not the happiest time complexity

    SUMMARY:
    I'm not quite sure yet
    We'll see what happens

'''
import urllib2 


def LZ78(s):    				# standard LZ78 algorithm
    s=s.upper()     				# make the sequence uppercase
    D={'A':0,'C':0,'G':0,'T':0}    		# initial dictionary
    cursor_pos=0    				# place the cursor to the beginning

    while cursor_pos < len(s):        				# do for the whole sequence
        word_len=1                  				# start with 1 letter word
        while s[cursor_pos:cursor_pos+word_len] in D: 		#if we have this word in the dictionary
            if (cursor_pos+word_len)==len(s):    		#if we reached the end stop the process
                break
            D[s[cursor_pos:cursor_pos+word_len]] +=1  		# counting how many times we look for this word in the dictionary
            word_len +=1            				# we had the word in the dictionary, now look for a longer one with an extra letter at the end
        if (cursor_pos+word_len)<=(len(s)):  			# if we are at the end of the sequence, don't to anything
	    if (s[cursor_pos+word_len-1] != 'N') :
        	D[s[cursor_pos:cursor_pos+word_len]] =1   		# this is a new word, add it to the dictionarym and say it we encountered it once
        cursor_pos +=word_len   				#we added an item to the dictionary, now move to the next phrase, by changing the cursor position
    
    return D    				# report the dictionary to the program

def fna_read(name):
    fid = open(name)
    header = fid.readline()
    dna = fid.read().replace("\n","")
    dna = dna.replace("\r","")
    fid.close()
    return header, dna

def multifna_read(name):
    # TODO:
    #	Sanatize Inputs, all upper, if x is not(A or T or G or C), x = N
    fragList = []
    headerList = []
    fid = open(name)
    a = fid.readline()
    headerList.append(a)
    s = ''
    for line in fid:
        if line[0]== '>':
	    headerList.append(line.strip())
            fragList.append(s)
            s=''
        else:
            s = s+line.strip()
    fid.close()
    fragList.append(s)
    return headerList, fragList

def main_getDEWList(classA, classB):
    '''
    OK, so first, I need to get each sequence, and for each sequence, generate a dictionary.
    
    First output goal:
	.ods with all words in each row,
	and all the columns with the specific sequence's counts of those words

    '''
    outfileName		= "outfile.ods"
    masterDict 		= {}
    masterWords 	= []
    classANames 	= []
    classASeqs		= []
    classADicts 	= []
    classBNames 	= []
    classBSeqs		= []
    classBDicts 	= []

    classANames, classASeqs = multifna_read(classA)
    classBNames, classBSeqs = multifna_read(classB)

    for i in range(len(classASeqs)):
	Di = LZ78(classASeqs[i])
	classADicts.append(Di)
	for word in Di.keys():
	    if word not in masterDict:
		masterDict[word] = 1
		masterWords.append(word)
    for j in range(len(classBSeqs)):
	Dj = LZ78(classBSeqs[j])
	classBDicts.append(Dj)
	for word in Dj.keys():
	    if word not in masterDict:
		masterDict[word] = 1
		masterWords.append(word)
	
    
    print len(masterWords)
    print masterDict["A"]

    outfile = open(outfileName, 'w')
    outfile.write("Word\t")
    for i in range(len(classANames)):
	outfile.write(classANames[i].strip())
	outfile.write("\t")
    for j in range(len(classBNames)):
	outfile.write(classBNames[j].strip())
	outfile.write("\t")
    outfile.write("\n")

    for k in range(0, 1000):
	outfile.write(masterWords[k])
        outfile.write("\t")
	for i in range(len(classADicts)):
	    Di = classADicts[i]
	    if masterWords[k] in Di:
		outfile.write(str(Di[masterWords[k]]))
	    else:
		outfile.write("0")
	    outfile.write("\t")
	for j in range(len(classBDicts)):
	    Dj = classBDicts[j]
	    if masterWords[k] in Dj:
		outfile.write(str(Dj[masterWords[k]]))
	    else:
		outfile.write("0")
	    outfile.write("\t")
	outfile.write("\n")


    return

def parseGeneFile(geneFile):
    geneStart = []
    geneStop = []
    NCID = []
    gf = open(geneFile, 'r')
    for line in gf:
	if (line[0] == '>'):
	    # Parse this gene
	    # >lcl|NC_008724.1_gene_1 [gene=z001L] [locus_tag=ATCV1_z001L] [location=complement(99..344)]
	    # >lcl|NC_008724.1_gene_2 [gene=Z002R] [locus_tag=ATCV1_Z002R] [location=551..1327]
	    # >lcl|NC_008724.1_gene_53 [gene=Z053L] [locus_tag=ATCV1_Z053L] [location=complement(13838..14563)]
	    temp = line
	    temp = temp.split()
	    geneID = temp[1]
	    geneID = geneID.split("=")[1]
	    NCID.append(geneID[:-1])
	    loc = temp[len(temp)-1]
	    loc = loc.split("=")
	    loc = loc[1]
	    if (loc[0] == 'c'):
		loc = loc[11:-2]
	    else:
		loc = loc[:-1]
	    loc = loc.split("..")
	    geneStart.append(int(loc[0]))
	    geneStop.append(int(loc[1]))
	    

    return geneStart, geneStop, NCID

def checkNCBI(geneID):
    target_url = "http://www.ncbi.nlm.nih.gov/gene/?term=" + str(geneID) + "&report=docsum&format=text"
    geneInfo = []
    for line in urllib2.urlopen(target_url):
	if (line[0] != "<"):
	    #print line.strip()
	    geneInfo.append(line.strip())
    #print len(geneInfo)
    #print geneInfo
    return geneInfo

def regionFoundIn(word, seq, geneFile):

    geneStart, geneStop, geneID = parseGeneFile(geneFile)
    #print len(geneStart), len(geneStop), len(geneID)
    name = seq[0].strip()
    seq = seq[1]
    word = word.upper()
    hits 		= 0
    geneRegion 		= 0
    grFlag		= False
    intergenicRegion 	= 0
    straddled 		= 0
    geneInfo 		= []
    genes		= []
    tempGenes		= []
    geneOut		= 'geneOut.txt'
    uniqGene		= {}
    uniqGeneList	= []
    uniqGeneID		= []

    for i in range(len(seq) - len(word)):
	checkWord = seq[i:i+len(word)].upper()
	if (checkWord == word):
	    grFlag = False
	    
	    for j in range(len(geneStart)):
		if ( (geneStart[j] <= i) and (i+len(word) <= geneStop[j]) ):
		    grFlag = True
		    geneRegion += 1
		    hits += 1
		    #print geneID[j]
		    try:
			geneInfo.append(checkNCBI(geneID[j]))
		    except:
			print "Opps! Bad URL Path!\n\tChech the Tubes!"
			
	    if (grFlag == False):
		intergenicRegion += 1
		hits += 1
	
	
    for i in range(len(geneInfo)):
	j = 0
	#print len(geneInfo[i])
	#print geneInfo[i]
	while True:
	#for j in range(len(geneInfo[i])):
	    if (j >= len(geneInfo[i])-2):
		break
	    #print j
	    #print geneInfo[i][j+1]
	    a = geneInfo[i][j+1].split()
	    if (a[0] != "hypothetical"):
		genes.append(geneInfo[i][j+1])
		temp = geneInfo[i][j].split()[1]
		tempGenes.append(geneInfo[i][j].split()[1])
	    j += 6
	    
 
    print "Word:\t\t\t\t\t", word
    print "Total hits:\t\t\t\t", hits
    print "Total in gene region:\t\t\t", geneRegion
    print "Total in intergenic regions:\t\t", intergenicRegion
    if (hits == 0):
	hits = 1
    print "Percent of time in gene region:\t\t", 	geneRegion/float(hits)
    print "Percent of time in intergenic regions:\t", intergenicRegion/float(hits)
    print "Non Hypothetical Genes:\t\t\t", len(genes)	
    for i in range(len(genes)):
	print "\t\t\t\t\t", tempGenes[i], "\t", genes[i] 
    
		
    outfile = open(geneOut, 'a')
    outfile.write(str(name))
    outfile.write("\nWord:\t\t\t\t\t" + str(word))
    outfile.write("\nTotal hits:\t\t\t\t" + str(hits))
    outfile.write("\nTotal in gene region:\t\t\t" + str(geneRegion))
    outfile.write("\nTotal in intergenic regions:\t\t" + str(intergenicRegion))
    outfile.write("\nPercent of time in gene region:\t\t" + str(geneRegion/float(hits)))
    outfile.write("\nPercent of time in intergenic regions:\t" + str(intergenicRegion/float(hits)))
    for i in range(len(genes)):
	if genes[i] not in uniqGene:
	    uniqGene[genes[i]] = 1
	    uniqGeneList.append(genes[i])
	    uniqGeneID.append(tempGenes[i])

    outfile.write("\nNon Hypothetical Genes:\t\t\t" + str(len(uniqGeneList)))
    for i in range(len(uniqGeneList)):
	outfile.write("\n\t\t\t\t\t" + str(uniqGeneID[i]) + "\t" + str(uniqGeneList[i]))
	
    outfile.write("\n")
    outfile.close()

    '''
    Word
    Percent of time in gene region
    Percent of time in intergenic regions
    Percent of time in both
    total hits    
    '''

def forArrayofWords(word):

    seqF = "Acanthocystis_turfacea_Chlorella_virus_1_complete_genome.tfa"
    geneFile = "Acanthocystis_turfacea_Chlorella_GeneList.txt"
    seq = fna_read(seqF)
    print "Sequence: ",seq[0].strip()
    regionFoundIn(word, seq, geneFile)


    seqF = "Paramecium_bursaria_Chlorella_virus_1_complete_genome.tfa"
    geneFile = "Paramecium_bursaria_Chlorella_GeneList.txt"
    seq = fna_read(seqF)
    print "Sequence: ",seq[0].strip()
    regionFoundIn(word, seq, geneFile)
    print "\n"

#*********************************************************************#
classAMulti = "Acanthocystis_turfacea_Chlorella.tfa"
classBMulti = "Paramecium_bursaria_Chlorella.tfa"
#main_getDEWList(classAMulti, classBMulti)
#The above generates the words, manipulate in excel

outfile = open('geneOut.txt', 'w')
outfile.close()
words = ['GCTTCCCG','TTGTGATT','CTCCTCGT','CGTCCTCC','CAATTAT','CCGGGAGA','CATACGCT','CGTCCTCG','CAGAGGCC','ATTGTTTC','GGAGCAC','TAATTGTA','CTCCTCGG','ATCATAG','ATTGTTTG','TTGTTGAAC','TGGCAGGGG','CCTGGAC','ATCATAA','CATACGCG','CAGAGGCG','CGGCACG','CAGAGGCT']

#words = ["TTGTGATT", "TCGTAGCCG", "CCTCCGTG", "CATGGAGG", "TTGTGATT", "GCTTCCCG", "GAGACAT", "ATGTTA", "ATCATAA", "TCACAA", "AACATT", "CACGGT" ]
for i in range(len(words)):
    forArrayofWords(words[i])
    outfile = open('geneOut.txt', 'a')
    outfile.write("\n")
    outfile.close()






















