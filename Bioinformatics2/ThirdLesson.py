# genCode = {"A": ("GCA", "GCU", "GCG", "GCC"), "G": ("GGA", "GGU", "GGG", "GGC"), "V": ("GUA", "GUU", "GUC", "GUG"),
##           "E": ("GAA", "GAG"), "D": ("GAU", "GAC"),
#           "P": ("CCA", "CCU", "CCG", "CCC"), "R": ("CGA", "CGU", "CGG", "CGC", "AGG", "AGA"),
#           "L": ("CUA", "CUU", "CUC", "CUG", "UUG", "UUA"), "Q": ("CAA", "CAG"), "H": ("CAU", "CAC"),
#           "S": ("UCA", "UCU", "UCG", "UCC", "AGC", "AGU"), "T": ("ACA", "ACU", "ACG", "ACC"), "F": ("UUU", "UUC"),
#           "W": ("UGG"),
#           "C": ("UGU", "UGC"), "Y": ("UAC", "UAU"), "I": ("AUA", "AUU", "AUC"), "M": ("AUG",),
#           "K": ("AAA", "AAG"), "N": ("AAC", "AAU")
#           }
import itertools


def translateRNA(RNAseq):
    basePerCodon = 3
    genCode = {"AAA": "K", "AAU": "N", "AAG": "K", "AAC": "N",
               "AUA": "I", "AUU": "I", "AUC": "I", "AUG": "M",
               "AGA": "R", "AGU": "S", "AGG": "R", "AGC": "S",
               "ACA": "T", "ACU": "T", "ACG": "T", "ACC": "T",

               "UAA": "stop", "UAU": "Y", "UAG": "stop", "UAC": "Y",
               "UUA": "L", "UUU": "F", "UUC": "F", "UUG": "L",
               "UGA": "stop", "UGU": "C", "UGG": "W", "UGC": "C",
               "UCA": "S", "UCU": "S", "UCG": "S", "UCC": "S",

               "CAA": "Q", "CAU": "H", "CAG": "Q", "CAC": "H",
               "CUA": "L", "CUU": "L", "CUC": "L", "CUG": "L",
               "CGA": "R", "CGU": "R", "CGG": "R", "CGC": "R",
               "CCA": "P", "CCU": "P", "CCG": "P", "CCC": "P",

               "GAA": "E", "GAU": "D", "GAG": "E", "GAC": "D",
               "GUA": "V", "GUU": "V", "GUC": "V", "GUG": "V",
               "GGA": "G", "GGU": "G", "GGG": "G", "GGC": "G",
               "GCA": "A", "GCU": "A", "GCG": "A", "GCC": "A",
               }
    protSeq = list()
    for i in range(0, len(RNAseq) - basePerCodon + 1, 3):
        codon = RNAseq[i: i + basePerCodon]
        if genCode[codon] == "stop":
            return "".join(protSeq)
        protSeq.append(genCode[codon])
    return "".join(protSeq)


def reverseComplement(text):
    DNAdict = {"A": "T", "C": "G", "T": "A", "G": "C"}
    res = ""
    for letter in text:
        res += DNAdict[letter]
    return res[-1::-1]


# Написать генератор для получения кодонов

def findDNASubstringCodingGivenProt(protSeq, DNAseq):
    basePerCodon = 3
    substrings = []
    genCodeDNA = {"AAA": "K", "AAT": "N", "AAG": "K", "AAC": "N",
                  "ATA": "I", "ATT": "I", "ATC": "I", "ATG": "M",
                  "AGA": "R", "AGT": "S", "AGG": "R", "AGC": "S",
                  "ACA": "T", "ACT": "T", "ACG": "T", "ACC": "T",

                  "TAA": "stop", "TAT": "Y", "TAG": "stop", "TAC": "Y",
                  "TTA": "L", "TTT": "F", "TTC": "F", "TTG": "L",
                  "TGA": "stop", "TGT": "C", "TGG": "W", "TGC": "C",
                  "TCA": "S", "TCT": "S", "TCG": "S", "TCC": "S",

                  "CAA": "Q", "CAT": "H", "CAG": "Q", "CAC": "H",
                  "CTA": "L", "CTT": "L", "CTC": "L", "CTG": "L",
                  "CGA": "R", "CGT": "R", "CGG": "R", "CGC": "R",
                  "CCA": "P", "CCT": "P", "CCG": "P", "CCC": "P",

                  "GAA": "E", "GAT": "D", "GAG": "E", "GAC": "D",
                  "GTA": "V", "GTT": "V", "GTC": "V", "GTG": "V",
                  "GGA": "G", "GGT": "G", "GGG": "G", "GGC": "G",
                  "GCA": "A", "GCT": "A", "GCG": "A", "GCC": "A",
                  }
    revDNAseq = reverseComplement(DNAseq)
    for i in range(len(DNAseq) - len(protSeq) * basePerCodon):
        codon = DNAseq[i: i + basePerCodon]
        revCodon = revDNAseq[i: i + basePerCodon]
        if genCodeDNA[codon] == protSeq[0]:
            j = i
            for j in range(i, i + len(protSeq) * basePerCodon, basePerCodon):
                codon = DNAseq[j: j + basePerCodon]
                if genCodeDNA[codon] != protSeq[(j - i) // basePerCodon]:
                    break
            else:
                substrings.append(DNAseq[i: i + len(protSeq) * basePerCodon])
        if genCodeDNA[revCodon] == protSeq[0]:
            j = i
            for j in range(i, i + len(protSeq) * basePerCodon, basePerCodon):
                codon = revDNAseq[j: j + basePerCodon]
                if genCodeDNA[codon] != protSeq[(j - i) // basePerCodon]:
                    break
            else:
                substrings.append(DNAseq[-(i + len(protSeq) * basePerCodon): -i])
    return substrings


def calcProtMass(protSeq):
    aminoacidMass = {"G": 57, "A": 71, "S": 87, "P": 97, "V": 99, "T": 101, "C": 103, "I": 113, "L": 113, "N": 114,
                     "D": 115, "K": 128, "Q": 128, "E": 129, "M": 131, "H": 137, "F": 147, "R": 156, "Y": 163, "W": 186}
    mass = 0
    for aminoacid in protSeq:
        mass += aminoacidMass[aminoacid]
    return mass

def calcProtMass_v3(protSeq): #'127 134 156'
    return sum(list(map(int, protSeq.split())))

def massSpecCircularFragmentsMass(protSeq):  # Исправить название
    fragments = {protSeq: 1, "": 1}
    doubledProtSeq = protSeq * 2
    for fragmentLen in range(1, len(protSeq)):
        for j in range(len(protSeq)):
            fragment = doubledProtSeq[j: j + fragmentLen]
            if fragment not in fragments:
                fragments[fragment] = 0
            fragments[fragment] += 1
    return {calcProtMass(fragment): fragments[fragment] for fragment in fragments}


def massSpecLinearFragmentsMass(protSeq):  # Исправить название
    fragments = {protSeq: 1, "": 1}
    for fragmentLen in range(1, len(protSeq)):
        for j in range(len(protSeq) - fragmentLen + 1):
            fragment = protSeq[j: j + fragmentLen]
            if fragment not in fragments:
                fragments[fragment] = 0
            fragments[fragment] += 1
    return {calcProtMass(fragment): fragments[fragment] for fragment in fragments}


def genSpectrumForCircularPeptide(protSeq):
    spectrum = {str(calcProtMass(protSeq)): 1, '0': 1}
    doubledProtSeq = protSeq * 2
    for fragmentLen in range(1, len(protSeq)):
        for j in range(len(protSeq)):
            fragment = doubledProtSeq[j: j + fragmentLen]
            fragmentMass = str(calcProtMass(fragment))
            if fragmentMass not in spectrum:
                spectrum[fragmentMass] = 0
            spectrum[fragmentMass] += 1
    return spectrum

def genSpectrumForCircularPeptide2(protSeq):
    spectrum = {calcProtMass(protSeq): 1, 0: 1}
    doubledProtSeq = protSeq * 2
    for fragmentLen in range(1, len(protSeq)):
        for j in range(len(protSeq)):
            fragment = doubledProtSeq[j: j + fragmentLen]
            fragmentMass = calcProtMass(fragment)
            if fragmentMass not in spectrum:
                spectrum[fragmentMass] = 0
            spectrum[fragmentMass] += 1
    return spectrum

def genSpectrumForCircularPeptide_v3(protSeq):
    protSeq = (list(map(int, protSeq.split())))
    spectrum = {sum(protSeq): 1, 0: 1}
    protLen = len(protSeq)
    for fragmentLen in range(1, protLen):
        for j in range(protLen):
            if j + fragmentLen < protLen:
                fragment = protSeq[j: j + fragmentLen]
            else:
                fragment = protSeq[j:] + protSeq[:j + fragmentLen - protLen]
            fragmentMass = sum(fragment)
            if fragmentMass not in spectrum:
                spectrum[fragmentMass] = 0
            spectrum[fragmentMass] += 1
    return spectrum

def genSpectrumForLinearPeptide(protSeq):
    spectrum = {calcProtMass(protSeq): 1, 0: 1}
    for fragmentLen in range(1, len(protSeq)):
        for j in range(len(protSeq) - fragmentLen + 1):
            fragment = protSeq[j: j + fragmentLen]
            fragmentMass = calcProtMass(fragment)
            if fragmentMass not in spectrum:
                spectrum[fragmentMass] = 0
            spectrum[fragmentMass] += 1
    return spectrum

def genSpectrumForLinearPeptide_v3(protSeq):
    protSeq = (list(map(int, protSeq.split())))
    spectrum = {sum(protSeq): 1, 0: 1}
    protLen = len(protSeq)
    for fragmentLen in range(1, protLen):
        for j in range(protLen - fragmentLen + 1):
            fragment = protSeq[j: j + fragmentLen]
            fragmentMass = sum(fragment)
            if fragmentMass not in spectrum:
                spectrum[fragmentMass] = 0
            spectrum[fragmentMass] += 1
    return spectrum
print('________')#"D": 115, "K": 128, "Q": 128, "E": 129, "M": 131, "H": 137
print(genSpectrumForCircularPeptide_v3('115 128 131'))
print(genSpectrumForLinearPeptide_v3('128 128 128'))
print(genSpectrumForCircularPeptide('DQM'))
print(genSpectrumForLinearPeptide('KKK'))
print('________')

def findPeptideGivenMass(mass):
    aminoacids = "GASPVTCILNDKQEMHFRYW"


# returns list of all peptides that is longer then the given peptide by one aminoacid
def expandPeptide(peptide):
    aminoacids = "GASPVTCILNDKQEMHFRYW"
    return [peptide + aminoacid for aminoacid in aminoacids]


def expandPeptides(peptides):
    aminoacids = "GASPVTCILNDKQEMHFRYW"
    expandedPeptides = set()
    if len(peptides) == 0:
        return set(aminoacids)
    for peptide in peptides:
        expandedPeptides |= {peptide + aminoacid for aminoacid in aminoacids}
    return expandedPeptides

def expandPeptides_v3(peptides):
    masses = {'57', '71', '87', '97', '99', '101', '103', '113', '114', '115', '128', '129', '131', '137', '147', '156', '163', '186'}
    expandedPeptides = set()
    if len(peptides) == 0:
        return masses
    for peptide in peptides:
        expandedPeptides |= {peptide + ' ' + aminoacid for aminoacid in masses}
    return expandedPeptides

def expandPeptides_v4(peptides, aminoacids):
    expandedPeptides = set()
    if len(peptides) == 0:
        return set(aminoacids)
    for peptide in peptides:
        expandedPeptides |= {peptide + ' ' + aminoacid for aminoacid in aminoacids}
    return expandedPeptides

print(expandPeptides_v3({'127', '128'}))



def isSubDict(potentialSubDict, dict):  # set doesn't fit because of repeating values
    for k in potentialSubDict:
        if k not in dict:
            return False
        elif potentialSubDict[k] > dict[k]:
            return False
    return True


def isEqualDict(potentialEqualDict, dict):
    for k in potentialEqualDict:
        if k not in dict:
            return False
        elif potentialEqualDict[k] != dict[k]:
            return False
    return True


def makeDictFromSpectrum(spectrum):
    spectrumDict = dict()
    for fragment in spectrum:
        if fragment not in spectrumDict:
            spectrumDict[fragment] = 0
        spectrumDict[fragment] += 1
    return spectrumDict


def replaceAminoacidsWMasses(protSeq):
    aminoacidMass = {"G": 57, "A": 71, "S": 87, "P": 97, "V": 99, "T": 101, "C": 103, "I": 113, "L": 113, "N": 114,
                     "D": 115, "K": 128, "Q": 128, "E": 129, "M": 131, "H": 137, "F": 147, "R": 156, "Y": 163, "W": 186}
    result = "-".join([str(aminoacidMass[aminoacid]) for aminoacid in protSeq]).lstrip()
    return result


# Branch&Bound algorithm
def findPeptidesThatGiveConsistentSpectrum(
        spectrum):  # Звменить добавление аминокислот массами, чтобы не гонять разные наборы аминокислот с одной массой
    # Так же можно смотреть из каких амрнокислот состоит  белок, так как масс спектр идеальный и перебирать только их
    spectrumDict = makeDictFromSpectrum(spectrum)
    potentialPeptides = set("")
    notConsistentPeptides = set()
    peptides = set()
    maxMass = max(*spectrum)
    while True:
        potentialPeptides = expandPeptides(potentialPeptides)
        for peptide in potentialPeptides:
            if calcProtMass(peptide) == maxMass:
                if isEqualDict(genSpectrumForCircularPeptide(peptide), spectrumDict):
                    # if genSpectrumForCircularPeptide(peptide) == spectrumDict:
                    peptides.add(peptide)
                    notConsistentPeptides.add(peptide)
            elif not isSubDict(genSpectrumForLinearPeptide(peptide), spectrumDict):
                notConsistentPeptides.add(peptide)
        potentialPeptides -= notConsistentPeptides
        if len(potentialPeptides) == 0:
            break
    return {replaceAminoacidsWMasses(peptide) for peptide in peptides}


def findPeptidesThatGiveConsistentSpectrum2(spectrum):  # Звменить добавление аминокислот массами, чтобы не гонять разные наборы аминокислот с одной массой
    # Так же можно смотреть из каких амрнокислот состоит  белок, так как масс спектр идеальный и перебирать только их
    spectrumDict = makeDictFromSpectrum(spectrum)
    potentialPeptides = set("")
    notConsistentPeptides = set()
    peptides = set()
    maxMass = max(*spectrum)
    while True:
        potentialPeptides = expandPeptides(potentialPeptides)
        for peptide in potentialPeptides:
            if calcProtMass(peptide) == maxMass:
                if isEqualDict(genSpectrumForCircularPeptide(peptide), spectrumDict):
                    # if genSpectrumForCircularPeptide(peptide) == spectrumDict:
                    peptides.add(peptide)
                    notConsistentPeptides.add(peptide)
            elif not isSubDict(genSpectrumForLinearPeptide(peptide), spectrumDict):
                notConsistentPeptides.add(peptide)
        potentialPeptides -= notConsistentPeptides
        if len(potentialPeptides) == 0:
            break
    return peptides
