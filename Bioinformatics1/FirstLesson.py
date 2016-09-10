import urllib.request as url, timeit, itertools

def findPattern(text, pat):
    res = text.find(pat)
    if res == -1:
        return 0
    else:
        return 1 + findPattern(text[res + 1:], pat)
print(findPattern("CGCGATACGTTACATACATGATAGACCGCGCGCGATCATATCGCGATTATC", "CGCG"))

def findMostFreqkmers(text, k):
    kmers = dict()
    maxCount = 0
    mostFreqKmers = set()
    for i in range(len(text) - k):
        kmer = text[i:i+k]
        #(kmer)
        if kmer not in kmers:
            kmers[kmer] = findPattern(text[i:], kmer)
        if kmers[kmer] > maxCount:
            maxCount = kmers[kmer]
            mostFreqKmers.clear()
            mostFreqKmers.add(kmer)
            #print(mostFreqKmers)
        elif kmers[kmer] == maxCount:
            mostFreqKmers.add(kmer)
    return mostFreqKmers

print(findMostFreqkmers("CGCCTAAATAGCCTCGCGGAGCCTTATGTCATACTCGTCCT", 3))

def findkmersWithKnownFreq(text, k, t):
    kmers = dict()
    for i in range(len(text) - k):
        kmer = text[i:i+k]
        if kmer not in kmers:
            kmers[kmer] = findPattern(text[i:], kmer)
            if kmers[kmer] != t:
                kmers.pop(kmer)
    return kmers

def reverseComplement(text):
    DNAdict = {"A":"T","C":"G","T":"A","G":"C"}
    res = ""
    for letter in text:
        res += DNAdict[letter]
    return res[-1::-1]

print(reverseComplement("GCTAGCT"))
# CODE CHALLENGE: Solve the Pattern Matching Problem.
#    Input: Two strings, Pattern and Genome.
#    Output: A collection of space - separated integers specifying all starting positions where Pattern appears
#    as a substring of Genome.

def findPatternStartPositions(pat, text):
    startPos = []
    shift = 0
    ind = text.find(pat)
    while ind != -1:
        startPos.append(ind)
        shift = startPos[-1] + 1
        ind = text.find(pat, shift)
    return startPos

#with url.urlopen("https://stepic.org/media/attachments/lessons/3/Vibrio_cholerae.txt") as file:
#    data = str(file.read().strip())
#with open("C:/Users/Alex1_000/Desktop/Vibrio_cholerae.txt", "r") as data:
#    data = data.read().strip()

# Clump Finding Problem: Find patterns forming clumps in a string.
#   Input: A string Genome, and integers k, L, and t. Output: All distinct k - mers forming(L, t) - clumps in Genome.
def findClumpOfkmer(text, k, L, t):
    kmersdict = {}
    for i in range(len(text) - L + 1):  #Is 1 needed?
        part = text[i:i+L]
        kmersdict.update(findkmersWithKnownFreq(part, k, t))
        print(kmersdict)

    return kmersdict

def patternToNumber(pat):
    total = 0
    nucleotide_to_number = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i in range(len(pat)):
        total = total * 4 + nucleotide_to_number[pat[i]]
    return total

print(patternToNumber("GCTCTTGCATACAACCATA"))

def numberToPattern(num, k):
    pat = []
    number_to_nucleotide = {0:'A', 1:'C', 2:'G', 3:'T'}
    while k != 0:
        letter = number_to_nucleotide[num%4]
        pat.append(letter)
        num //= 4
        k -= 1
    #pat.append(number_to_nucleotide[num])
    pat.reverse()
    return "".join(pat)

print(numberToPattern(5818, 11))

def frequencyArray(text, k):
    freqArray = [0 for i in range(4 ** k)]
    for i in range(len(text) - k + 1):
        pat = text[i:i+k]
        #print(pat)
        num = patternToNumber(pat)
        freqArray[num] += 1
    return freqArray


def findingFrequentWordsBySorting(text , k):
    freqSet = ()
    index = []
    count = []
    for i in range(len(text) - k + 1):
        pat = text[i:i + k]
        num = patternToNumber(pat)
        index.append(patternToNumber(pat))
        count.append(1)
    sortedIndex = index.sort()
    print(sortedIndex)
    for i in range(1, len(text) - k + 1):
        if sortedIndex[i] == sortedIndex[i - 1]:
            count[i] = count[i - 1] + 1
    maxCount = max(*count)
    for i in range(len(text) - k + 1):
        if count[i] == maxCount:
            pat = numberToPattern(sortedIndex[i], k)
            freqSet.add(pat)

    return freqSet

#text = "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA"
#k, L, t = 5, 50, 4
#for i in findClumpOfkmer(text, k, L, t).keys():
#    print(i, end=" ")


def betterClumpFinding(genome, k, L, t):
    frequentPatterns = set()
    clump = [0 for i in range(4 ** k)]
    text = genome[:L]
    freqArray = frequencyArray(text, k)
    for i in range(4 ** k):
        if freqArray[i] >= t:
            clump[i] = 1
    for i in range(1, len(genome) - L + 1):
        firstPat = genome[i - 1: i - 1 + k]

        index = patternToNumber(firstPat)
        freqArray[index] -= 1
        lastPat = genome[(i + L - k): i + L]

        index = patternToNumber(lastPat)
        freqArray[index] += 1
        if freqArray[index] >= t:
            clump[index] = 1
    for i in range(4**k):
        if clump[i] == 1:
            pat = numberToPattern(i, k)
            frequentPatterns.add(pat)
    print(frequentPatterns)
    return frequentPatterns
text = "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA"
k, L, t = 5, 50, 4
for i in betterClumpFinding(text, k, L, t):
    print(i, end=" ")