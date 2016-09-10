# We define Skewi(Genome) as the difference between the total number of occurrences of G and the total number of occurrences
#  of C in the first i nucleotides of Genome. The skew diagram is defined by plotting Skewi (Genome) (as i ranges from 0 to
# |Genome|), where Skew0 (Genome) is set equal to zero.
def skew(seq, n):
    skewAr = [0]
    f = lambda l: -1 if (l == "C") else 1 if (l == "G") else 0
    for i in seq[:n]:
        skewAr.append(f(i) + skewAr[-1])
    return skewAr


# print(*skew("GAGCCACCGCGATA", 14))

# Minimum Skew Problem: Find a position in a genome where the skew diagram attains a minimum
# Input: A DNA string Genome
# Output: All integer(s) i minimizing Skewi (Genome) among all values of i (from 0 to |Genome|)
def minimumSkew(seq):
    skewAr = skew(seq, len(seq))
    print(skewAr)
    min = 0
    minSkewAr = []
    for i in range(len(skewAr)):
        if skewAr[i] < min:
            min = skewAr[i]
            minSkewAr = [i]
        elif skewAr[i] == min:
            minSkewAr.append(i)
    return minSkewAr

print(*minimumSkew("CATTCCAGTACTTCGATGATGGCGTGAAGA"))

#The number of mismatches between strings p and q is called the Hamming distance between these strings and
# is denoted HammingDistance(p, q)
#Hamming Distance Problem: Compute the Hamming distance between two strings
#Input: Two strings of equal length
#Output: The Hamming distance between these strings
def hammingDist(seq1, seq2):
    assert len(seq1) == len(seq2)
    hamDist = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            hamDist += 1
    return hamDist

print("Ham:", hammingDist("CTACAGCAATACGATCATATGCGGATCCGCAGTGGCCGGTAGACACACGT","CTACCCCGCTGCTCAATGACCGGGACTAAAGAGGCGAAGATTATGGTGTG"))

#Approximate Pattern Matching Problem: Find all approximate occurrences of a pattern in a string
#Input: Strings Pattern and Text along with an integer d
# Output: All starting positions where Pattern appears as a substring of Text with at most d mismatches
def findApproxMatchesSlow(pat, text, er):
    startPos = []
    ind = 0
    while ind <= (len(text) - len(pat)):
        match = True
        mism = 0
        for i in range(len(pat)):
            if mism > er:
                match = False
                break
            #print(pat[i], text[ind + i])
            if pat[i] != text[ind + i]:
                mism += 1
            if mism > er:
                match = False
                break
        if match:
            startPos.append(ind)
        ind += 1
        #print(ind)
    return startPos
#pat = "ATATATCTGA"

# er = 5
# print(*findApproxMatchesSlow(pat, text, er))


def countApproxPatternMathes(pat, text, er):
    return len(findApproxMatchesSlow(pat, text, er))

print(countApproxPatternMathes("TGT", "CGTGACAGTGTATGGGCATCTTT", 1))


def neighbours(pat, er):
    nucl = ("A", "C", "G", "T")
    if er == 0:
        neighborhood = set()
        neighborhood.add(pat)
        #print(neighborhood)
        return neighborhood
    if len(pat) == 1:
        return nucl
    neighborhood = set()
    suffixNeighbours = neighbours(pat[1:], er)
    for neighbor in suffixNeighbours:
        if hammingDist(pat[1:], neighbor) < er:
            for nucleotide in nucl:
                neighborhood.add(nucleotide + neighbor)
        else:
            neighborhood.add(pat[0] + neighbor)
    return neighborhood

print(*(neighbours("TAAAAAAAAGGGGGG", 1)))
#print(*neighbours("GGCCCAGAG", 3))
#for i in neighbours("GTACGGCC", 2):
#    print(i)


def findMostFreqkmersWMismatches(text, k, d):
    mostFreqKmers = set()
    kmersDict = {}
    maxCount = 0
    for i in range(len(text) - k):
        kmer = text[i:i + k]
        print(kmer)
        if kmer not in kmersDict:
            kmersDict[kmer] = countApproxPatternMathes(kmer, text, d)
            for neighbour in neighbours(kmer, d):
                #print(":", neighbour)
                kmersDict[kmer] += countApproxPatternMathes(neighbour, text, d)
        if kmersDict[kmer] > maxCount:
            maxCount = kmersDict[kmer]
            mostFreqKmers.clear()
            mostFreqKmers.add(kmer)
        elif kmersDict[kmer] == maxCount:
            mostFreqKmers.add(kmer)
    print(kmersDict)
    return mostFreqKmers



def findMostFreqkmersWMismatches2(text, k, d):
    # Use only with d > 0, otherwise neighboursDict will contain only one kmer and "for neighbour in neighbours"
    # will give  separated bases instead of kmers
    assert d > 0
    mostFreqKmers = set()
    kmersSet = set()
    neighboursDict ={}
    maxCount = 0
    text = text.upper()
    for i in range(len(text) - k):
        kmer = text[i:i + k]
        kmersSet.add(kmer)
        for neighbour in neighbours(kmer, d):
            if neighbour not in neighboursDict:
                neighboursDict[neighbour] = 0
            neighboursDict[neighbour] += 1
    for kmer in kmersSet:
        if neighboursDict[kmer] > maxCount:
            maxCount = neighboursDict[kmer]
            mostFreqKmers.clear()
            mostFreqKmers.add(kmer)
        elif neighboursDict[kmer] == maxCount:
            mostFreqKmers.add(kmer)

    return mostFreqKmers

def reverseComplement(text):
    DNAdict = {"A":"T","C":"G","T":"A","G":"C"}
    res = ""
    for letter in text:
        res += DNAdict[letter]
    return res[-1::-1]

def findMostFreqkmersWMismatches3(text, k, d):
    mostFreqKmers = set()
    kmersSet = set()
    neighboursDict ={}
    maxCount = 0
    for i in range(len(text) - k):
        kmer = text[i:i + k]
        kmersSet.add(kmer)
        kmersSet.add(reverseComplement(kmer))
        for neighbour in neighbours(kmer, d):
            #print(":", neighbour)
            if neighbour not in neighboursDict:
                neighboursDict[neighbour] = 0
            neighboursDict[neighbour] += 1
        for neighbour in neighbours(reverseComplement(kmer), d):
            # print(":", neighbour)
            if neighbour not in neighboursDict:
                neighboursDict[neighbour] = 0
            neighboursDict[neighbour] += 1
    for kmer in kmersSet:
        #print(kmer,":",neighboursDict[kmer])
        if neighboursDict[kmer] > maxCount:
            maxCount = neighboursDict[kmer]
            mostFreqKmers.clear()
            mostFreqKmers.add(kmer)
        elif neighboursDict[kmer] == maxCount:
            mostFreqKmers.add(kmer)

    return mostFreqKmers

print(findMostFreqkmersWMismatches3("CTATTCGGTTTTTTTTTTTCGTCTATTCGGCTATTCGTTTTCGGCGTCGTTTTCGGTTTTTTCTACTATTTTTTCGGTTTTCTATTTCGTTTTCGTCGGCGGTTTCGTCGGCGGCTACGGTTTTTCGGCGTCTACGGTTTTTCGTCGGCGTCGGCGTCGTCGTCGTTTTCGTCTATTTTCTACGGCTACGGTTTCGTCGGTTTCGTTTTTT", 5, 2))