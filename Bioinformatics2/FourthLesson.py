import Bioinf2.FirstLesson as first, Bioinf2.SecondLesson as sec, Bioinf2.ThirdLesson as third

def score(peptide, spectrum):
    score = 0
    spectrum_dict = third.makeDictFromSpectrum(spectrum)
    peptide_spectrum_dict = third.genSpectrumForLinearPeptide(peptide)
    print(spectrum_dict)
    print(peptide_spectrum_dict)
    for frag in peptide_spectrum_dict:
        if frag in spectrum_dict:
            score += min(peptide_spectrum_dict[frag], spectrum_dict[frag])
    return score
print('___!!!____')
print(score('PEEP', list(map(int,'0 97 129 129 129 194 226 323 323 355 452'.split()))))

def score2(peptide, spectrumDict):
    score = 0
    peptide_spectrum_dict = third.genSpectrumForCircularPeptide2(peptide)
    for frag in peptide_spectrum_dict:
        if frag in spectrumDict:
            score += min(peptide_spectrum_dict[frag], spectrumDict[frag])
    return score

def score3(peptide, spectrumDict):
    score = 0
    peptide_spectrum_dict = third.genSpectrumForLinearPeptide(peptide)
    for frag in peptide_spectrum_dict:
        if frag in spectrumDict:
            score += min(peptide_spectrum_dict[frag], spectrumDict[frag])
    return score

def cyclicScore_v3(peptide, spectrumDict):
    score = 0
    peptide_spectrum_dict = third.genSpectrumForCircularPeptide_v3(peptide)
    for frag in peptide_spectrum_dict:
        if frag in spectrumDict:
            score += min(peptide_spectrum_dict[frag], spectrumDict[frag])
    return score

def linearScore_v3(peptide, spectrumDict):
    score = 0
    peptide_spectrum_dict = third.genSpectrumForLinearPeptide_v3(peptide)
    for frag in peptide_spectrum_dict:
        if frag in spectrumDict:
            score += min(peptide_spectrum_dict[frag], spectrumDict[frag])
    return score

def score4(peptide, spectrum):
    score = 0
    spectrum_dict = third.makeDictFromSpectrum(spectrum)
    peptide_spectrum_dict = third.genSpectrumForLinearPeptide(peptide)
    print(spectrum_dict)
    print(peptide_spectrum_dict)
    for frag in peptide_spectrum_dict:
        if frag in spectrum_dict:
            score += min(peptide_spectrum_dict[frag], spectrum_dict[frag])
    return score

def leaderBoardCycloPeptideSequencing(n, spectrum):
    spectrumDict = third.makeDictFromSpectrum(spectrum)
    potentialPeptides = set("")
    leaderBoard = dict()
    allPeptide = dict()
    maxMass = max(*spectrum)
    leaderPeptide = ""
    while True:
        potentialPeptides = third.expandPeptides(potentialPeptides)
        for peptide in potentialPeptides:
            if peptide not in allPeptide:
                allPeptide[peptide] = score3(peptide, spectrumDict)
            leaderBoard[peptide] = allPeptide[peptide]
            pepMass = third.calcProtMass(peptide)
            if pepMass == maxMass:
                if score2(peptide, spectrumDict) > score2(leaderPeptide, spectrumDict):
                    leaderPeptide = peptide
                    #print(peptide)
                    print(third.replaceAminoacidsWMasses(leaderPeptide))
                del leaderBoard[peptide]
            elif pepMass > maxMass:
                del leaderBoard[peptide]
        potentialPeptides = set(sorted(leaderBoard, key=leaderBoard.get, reverse=True)[:n])
        print(potentialPeptides)
        leaderBoard = dict()
        if len(potentialPeptides) == 0:
            break
    return third.replaceAminoacidsWMasses(leaderPeptide)

def leaderBoardCycloPeptideSequencing_v2(n, spectrum):
    spectrumDict = third.makeDictFromSpectrum(spectrum)
    potentialPeptides = set("")
    leaderBoard = dict()
    allPeptide = dict()
    maxMass = max(*spectrum)
    leaderPeptide = ""
    leaderScore = 0
    while True:
        potentialPeptides = third.expandPeptides(potentialPeptides)
        for peptide in potentialPeptides:
            if peptide not in allPeptide:
                allPeptide[peptide] = score3(peptide, spectrumDict)
            leaderBoard[peptide] = allPeptide[peptide]
            pepMass = third.calcProtMass(peptide)
            if pepMass == maxMass:
                score = score2(peptide, spectrumDict)
                if score == 83:
                    print(peptide)
                if score > leaderScore:
                    leaderPeptide = peptide
                    #print(peptide)
                    #print(third.replaceAminoacidsWMasses(leaderPeptide))
            elif pepMass > maxMass:
                del leaderBoard[peptide]
        potentialPeptides = sorted(leaderBoard, key=leaderBoard.get, reverse=True)
        j = n - 1
        while ((len(potentialPeptides) - 1) > j) and (leaderBoard[potentialPeptides[j + 1]] == leaderBoard[potentialPeptides[j]]): #since potentia;Peptides doesn't store score
            j += 1
        potentialPeptides = potentialPeptides[:j + 1]
        print('!')
        #163 - 147 - 129 - 137 - 128 - 115 - 128 - 103 - 156 - 97 - 163 - 99 - 103 - 128 - 163 - 163 - 99 - 103 - 99 - 186

        if len(potentialPeptides) == 0:
            #print(leaderBoard)
            break
        leaderBoard = dict()
    #for p in potentialPeptides:
    #    if allPeptide[p] == 38:
    #        print(p)
    return third.replaceAminoacidsWMasses(leaderPeptide), leaderPeptide

def trim(peptides, spectrum, n):
    leaderBoard = dict()
    spectrumDict = third.makeDictFromSpectrum(spectrum)
    for peptide in peptides:
        leaderBoard[peptide] = score3(peptide, spectrumDict)
    potentialPeptides = sorted(leaderBoard, key=leaderBoard.get, reverse=True)
    j = n - 1
    while ((len(potentialPeptides) - 2) > j) and (leaderBoard[potentialPeptides[j + 1]] == leaderBoard[potentialPeptides[j]]):  # since potentia;Peptides doesn't store score
        j += 1
    potentialPeptides = potentialPeptides[:j + 1 ]
    return potentialPeptides


def print_result(peptide):
    #print(*peptide)
    print('-'.join(peptide.split()).strip('-'))

def leaderBoardCycloPeptideSequencing_v3(n, spectrum):
    spectrumDict = third.makeDictFromSpectrum(spectrum)
    potentialPeptides = set("")
    leaderBoard = dict()
    allPeptide = dict()
    maxMass = max(*spectrum)
    leaderPeptide = ""
    leaderScore = 0
    while True:
        potentialPeptides = third.expandPeptides_v3(potentialPeptides)
        #print(potentialPeptides)
        for peptide in potentialPeptides:
            if peptide not in allPeptide:
                allPeptide[peptide] = linearScore_v3(peptide, spectrumDict)
            leaderBoard[peptide] = allPeptide[peptide]
            pepMass = third.calcProtMass_v3(peptide)
            if pepMass == maxMass:
                score = cyclicScore_v3(peptide, spectrumDict)
                if score == 83:
                    print('!!!')
                    print_result(peptide)
                if score > leaderScore:
                    leaderPeptide = peptide
            elif pepMass > maxMass:
                del leaderBoard[peptide]
        potentialPeptides = sorted(leaderBoard, key=leaderBoard.get, reverse=True)
        j = n - 1
        while ((len(potentialPeptides) - 1) > j) and (leaderBoard[potentialPeptides[j + 1]] == leaderBoard[potentialPeptides[j]]): #since potentia;Peptides doesn't store score
            j += 1
        potentialPeptides = potentialPeptides[:j + 1]
        print('!')
        if len(potentialPeptides) == 0:
            break
        leaderBoard = dict()

    return leaderPeptide



def convolution(spectrum):
    conv_spectrum = list()
    spectrum = sorted(spectrum)
    for i in range(len(spectrum) - 1, 0, -1):
        for j in range(i):
            delta = spectrum[i] - spectrum[j]
            if delta > 0:
                conv_spectrum.append(delta)
    return  conv_spectrum

def getPossibleAminoacidsByConvolution(m, spectrum):
    conv_spectrum_dict = dict()
    spectrum = sorted(spectrum)
    for i in range(len(spectrum) - 1, 0, -1):
        for j in range(i):
            delta = spectrum[i] - spectrum[j]
            if 57 <= delta <= 200:
                delta = str(delta)
                if delta not in conv_spectrum_dict:
                    conv_spectrum_dict[delta] = 0
                conv_spectrum_dict[delta] += 1
    conv_spectrum = sorted(conv_spectrum_dict, key=conv_spectrum_dict.get, reverse=True)
    j = m - 1
    while ((len(conv_spectrum) - 1) > j) and (conv_spectrum_dict[conv_spectrum[j + 1]] == conv_spectrum_dict[conv_spectrum[j]]):  # since potentia;Peptides doesn't store score
        j += 1
    conv_spectrum = conv_spectrum[:j + 1]
    return conv_spectrum

def convolutionCyclopeptideSequencing(m, n, spectrum):
    aminoacids = getPossibleAminoacidsByConvolution(m, spectrum)
    spectrumDict = third.makeDictFromSpectrum(spectrum)
    potentialPeptides = set("")
    leaderBoard = dict()
    allPeptide = dict()
    maxMass = max(*spectrum)
    leaderPeptide = ""
    leaderScore = 0
    max_score = 0
    while True:
        potentialPeptides = third.expandPeptides_v4(potentialPeptides, aminoacids)
        #print(potentialPeptides)
        for peptide in potentialPeptides:
            if peptide not in allPeptide:
                allPeptide[peptide] = linearScore_v3(peptide, spectrumDict)
            leaderBoard[peptide] = allPeptide[peptide]
            pepMass = third.calcProtMass_v3(peptide)
            if pepMass == maxMass:
                score = cyclicScore_v3(peptide, spectrumDict)
                max_score = max(score, max_score)
                if score == 82:
                    #print('!!!')
                    print_result(peptide)
                if score > leaderScore:
                    leaderPeptide = peptide
            elif pepMass > maxMass:
                del leaderBoard[peptide]
        potentialPeptides = sorted(leaderBoard, key=leaderBoard.get, reverse=True)
        j = n - 1
        while ((len(potentialPeptides) - 1) > j) and (leaderBoard[potentialPeptides[j + 1]] == leaderBoard[
            potentialPeptides[j]]):  # since potentia;Peptides doesn't store score
            j += 1
        potentialPeptides = potentialPeptides[:j + 1]
        print('!')
        if len(potentialPeptides) == 0:
            break
        leaderBoard = dict()
    print(max_score)
    #print_result(leaderPeptide)
    return leaderPeptide



print(convolution([0, 86, 160, 234, 308, 320, 382]))