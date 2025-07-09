# A collection of functions for DNA sequence analysis and motif finding.
# Credits: Bioinformatics For Beginners (UC San Diego, Coursera)
#
# This file contains implementations of various algorithms and helper functions
# for analyzing DNA sequences, including motif finding, pattern matching, and
# probabilistic models. Each function is documented with its purpose and usage.

import random

def PatternCount(Text, Pattern):
    """
    Counts the number of times a given pattern appears in a text.
    Args:
        Text (str): The DNA string to search within.
        Pattern (str): The pattern to search for.
    Returns:
        int: The number of occurrences of Pattern in Text.
    """
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 

def FrequentWords(Text, k):
    """
    Finds the most frequent k-mers (substrings of length k) in a DNA string.
    Args:
        Text (str): The DNA string to search within.
        k (int): The length of substrings to consider.
    Returns:
        list: A list of the most frequent k-mers in Text.
    """
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words

def FrequencyMap(Text, k):
    """
    Creates a frequency map of all k-mers in a DNA string.
    Args:
        Text (str): The DNA string to search within.
        k (int): The length of substrings to consider.
    Returns:
        dict: A dictionary mapping each k-mer to its frequency in Text.
    """
    freq = {}
    n = len(Text)
    
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
        for f in range(n-k+1):
            if Text[f:f+k] == Pattern:
                freq[Pattern] += 1
    return freq

def PatternMatching(Pattern, Genome):
    """
    Finds all starting positions where a pattern appears in a genome.
    Args:
        Pattern (str): The pattern to search for.
        Genome (str): The DNA string to search within.
    Returns:
        list: A list of starting positions where Pattern appears in Genome.
    """
    positions = []
    n = len(Pattern)
    
    for i in range(len(Genome) - n + 1):
        if Genome[i: i+n] == Pattern:
            positions.append(i)
    return positions

def ReverseComplement(Pattern):
    """
    Computes the reverse complement of a DNA string.
    Args:
        Pattern (str): The DNA string.
    Returns:
        str: The reverse complement of Pattern.
    """

    Pattern = Reverse(Pattern)
    Pattern = Complement(Pattern)
    return Pattern

def Reverse(Pattern):
    """
    Reverses a DNA string.
    Args:
        Pattern (str): The DNA string.
    Returns:
        str: The reversed DNA string.
    """
    return Pattern[::-1]  # Pythonic string reversal

def Complement(Pattern):
    """
    Computes the complement of a DNA string.
    Args:
        Pattern (str): The DNA string.
    Returns:
        str: The complement of Pattern.
    """
    compPattern = ""
    for char in Pattern:
        if char == "A": compPattern += "T"
        if char == "T": compPattern += "A"
        if char == "C": compPattern += "G"
        if char == "G": compPattern += "C"
    return compPattern

#####################################################################################################

def SymbolArray(Genome, symbol):
    """
    Computes the symbol array for a given symbol in a genome (inefficient O(n^2)).
    The symbol array is a list where each index i contains the count of the symbol
    in the first half of the genome starting from index i.
    Args:
        Genome (str): The DNA string to analyze.
        symbol (str): The symbol to count in the genome.
    Returns:
        dict: A dictionary where keys are indices and values are counts of the symbol.
    Note: This is a slow implementation for large genomes. See FasterSymbolArray.
    """
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array

def FasterSymbolArray(Genome, symbol):
    """
    Computes the symbol array for a given symbol in a genome efficiently (O(n)).
    Args:
        Genome (str): The DNA string to analyze.
        symbol (str): The symbol to count in the genome.
    Returns:
        dict: A dictionary where keys are indices and values are counts of the symbol.
    """
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        array[i] = array[i-1]
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

def MinimumSkew(Genome):
    """
    Finds all positions in the genome where the skew (difference between G and C)
    is minimized. Useful for locating the origin of replication.
    Args:
        Genome (str): The DNA string.
    Returns:
        list: Positions with minimum skew.
    """
    positions = []
    skew_list = SkewArray(Genome)
    min_skew = min(skew_list)
    for i,n in enumerate(skew_list):
        if n == min_skew:
            positions.append(i)
    return positions

def SkewArray(Genome):
    """
    Computes the skew array for a genome (difference between G and C at each position).
    Args:
        Genome (str): The DNA string.
    Returns:
        list: Skew values at each position.
    """
    skew = [0] 
    n = len(Genome)
    for i in range(n):
        if Genome[i] == "C":
            skew.append(skew[i] - 1)
        elif Genome[i] == "G":
            skew.append(skew[i] + 1)
        else:
            skew.append(skew[i])
    return skew

################################################################################################

def ApproximatePatternMatching(Text, Pattern, d):
    """
    Finds all positions where Pattern appears in Text with at most d mismatches.
    Args:
        Text (str): The DNA string to search within.
        Pattern (str): The pattern to search for.
        d (int): Maximum number of mismatches allowed.
    Returns:
        list: Starting positions with at most d mismatches.
    """
    positions = []
    n = len(Pattern)
    for i in range(len(Text)-n+1):
        if HammingDistance(Pattern, Text[i:i+n]) <= d:
            positions.append(i)
    return positions

def HammingDistance(p, q):
    """
    Computes the Hamming distance (number of mismatches) between two strings.
    Args:
        p (str): First string.
        q (str): Second string.
    Returns:
        int: Number of positions where p and q differ, or -1 if lengths differ.
    """
    len_p = len(p)
    hemm = 0
    if len_p != len(q):
        return -1
    else:
        for i in range(len_p):
            if p[i] != q[i]:
                hemm += 1
    return hemm

def ApproximatePatternCount(Pattern, Text, d):
    """
    Counts the number of times Pattern appears in Text with at most d mismatches.
    Args:
        Pattern (str): The pattern to search for.
        Text (str): The DNA string to search within.
        d (int): Maximum number of mismatches allowed.
    Returns:
        int: Number of approximate matches.
    """
    positions = []
    n = len(Pattern)
    for i in range(len(Text)-n+1):
        if HammingDistance(Pattern, Text[i:i+n]) <= d:
            positions.append(i)
    return len(positions)

###########################################################################

def Count_old(Motifs):
    """
    Counts the occurrences of each nucleotide at each position in a list of motifs (no pseudocounts).
    Args:
        Motifs (list): List of k-mer strings.
    Returns:
        dict: Nucleotide counts at each position.
    """
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
                
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
            
    return count

def Profile_old(Motifs):
    """
    Computes the profile matrix (probabilities) for a list of motifs (no pseudocounts).
    Args:
        Motifs (list): List of k-mer strings.
    Returns:
        dict: Profile matrix as probabilities.
    """
    t = len(Motifs)
    k = len(Motifs[0])
    profile = Count(Motifs)
    for i in 'ACGT':
        for j in range(k):
            profile[i][j] = profile[i][j]/t
    return profile

def Consensus(Motifs):
    """
    Computes the consensus string from a list of motifs (most common nucleotide at each position).
    Args:
        Motifs (list): List of k-mer strings.
    Returns:
        str: Consensus string.
    """
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
        
    return consensus

def Score(Motifs):
    """
    Computes the score of a set of motifs (total number of non-consensus nucleotides).
    Args:
        Motifs (list): List of k-mer strings.
    Returns:
        int: Score (lower is better).
    """
    consensus = Consensus(Motifs)
    k = len(consensus)
    t = len(Motifs)
    score = 0
    for j in range(k):
        for i in range(t):
            if Motifs[i][j] != consensus[j]:
                score += 1
    return score

def Pr(Text, Profile):
    """
    Computes the probability of a k-mer Text given a profile matrix.
    Args:
        Text (str): k-mer string.
        Profile (dict): Profile matrix.
    Returns:
        float: Probability of Text according to Profile.
    """
    p = 1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]
    return p

def ProfileMostProbableKmer(text, k, profile):
    """
    Finds the most probable k-mer in text according to a given profile matrix.
    Args:
        text (str): DNA string.
        k (int): Length of k-mer.
        profile (dict): Profile matrix.
    Returns:
        str: Most probable k-mer.
    """
    max_probability = -1
    most_probable_kmer = ""
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        probability = Pr(kmer, profile)
        if probability > max_probability:
            max_probability = probability
            most_probable_kmer = kmer
    return most_probable_kmer

def GreedyMotifSearch(Dna, k, t):
    """
    Implements the Greedy Motif Search algorithm for finding motifs in DNA sequences.
    Args:
        Dna (list): List of DNA strings.
        k (int): Length of motif.
        t (int): Number of DNA strings.
    Returns:
        list: Best motifs found.
    """
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

def Count(Motifs):
    """
    Counts the occurrences of each nucleotide at each position in a list of motifs (with pseudocounts).
    Args:
        Motifs (list): List of k-mer strings.
    Returns:
        dict: Nucleotide counts at each position (with pseudocounts).
    """
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)
                
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
            
    return count

def Profile(Motifs):
    """
    Computes the profile matrix (probabilities) for a list of motifs (with pseudocounts).
    Args:
        Motifs (list): List of k-mer strings.
    Returns:
        dict: Profile matrix as probabilities (with pseudocounts).
    """
    t = len(Motifs)
    k = len(Motifs[0])
    profile = Count(Motifs)
    for i in 'ACGT':
        for j in range(k):
            profile[i][j] = profile[i][j]/(t+4)
    return profile

def RandomMotifs(Dna, k, t):
    """
    Generates a random set of motifs (one from each DNA string).
    Args:
        Dna (list): List of DNA strings.
        k (int): Length of motif.
        t (int): Number of DNA strings.
    Returns:
        list: Randomly selected motifs.
    """
    random_motifs = []
    t = len(Dna)
    l = len(Dna[0])
    for i in range(t):
        ran_num = random.randint(0,l-k)
        random_motifs.append(Dna[i][ran_num:ran_num+k])
    return random_motifs

def RandomizedMotifSearch(Dna, k, t):
    """
    Implements the Randomized Motif Search algorithm.
    Args:
        Dna (list): List of DNA strings.
        k (int): Length of motif.
        t (int): Number of DNA strings.
    Returns:
        list: Best motifs found.
    """
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = Profile(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs 

def Motifs(pf,dna):
    """
    Finds the most probable motif in each DNA string given a profile matrix.
    Args:
        pf (dict): Profile matrix.
        dna (list): List of DNA strings.
    Returns:
        list: Most probable motifs for each string.
    """
    k = len(pf['A'])
    D = []
    for i in range(0,len(dna)):
        km = []
        sc = []
        for kk in range(len(dna[i])-k+1):
            km += [dna[i][kk:kk+k]]
        for i in km:
            sc += [Pr(i,pf)]
        D += [km[sc.index(max(sc))]]
    return D

def Normalize(Probabilities):
    """
    Normalizes a dictionary of probabilities so they sum to 1.
    Args:
        Probabilities (dict): Dictionary of k-mer probabilities.
    Returns:
        dict: Normalized probabilities.
    """
    total = sum(Probabilities.values())
    normalized = {kmer: prob / total for kmer, prob in Probabilities.items()}
    return normalized

def WeightedDie(Probabilities):
    """
    Selects a k-mer randomly according to the given probability distribution.
    Args:
        Probabilities (dict): Dictionary of k-mer probabilities.
    Returns:
        str: Randomly selected k-mer.
    """
    p = random.uniform(0, 1)
    for kmer, probability in Probabilities.items():
        p -= probability
        if p <= 0:
            return kmer

def ProfileGeneratedString(Text, profile, k):
    """
    Generates a k-mer from Text according to the profile matrix probabilities.
    Args:
        Text (str): DNA string.
        profile (dict): Profile matrix.
        k (int): Length of k-mer.
    Returns:
        str: Randomly generated k-mer.
    """
    n = len(Text)
    probabilities = {}
    for i in range(0, n - k + 1):
        kmer = Text[i:i+k]
        probabilities[kmer] = Pr(kmer, profile)
    probabilities = Normalize(probabilities) 
    return WeightedDie(probabilities)

def GibbsSampler(Dna, k, t, N):
    """
    Implements the Gibbs Sampler algorithm for motif finding.
    Args:
        Dna (list): List of DNA strings.
        k (int): Length of motif.
        t (int): Number of DNA strings.
        N (int): Number of iterations.
    Returns:
        list: Best motifs found.
    """
    Motifs=RandomMotifs(Dna, k, t)
    BestMotifs=Motifs
    for j in range(N):
        i=random.randint(1,t)
        text=Motifs.pop(i-1)
        profile=Profile(Motifs)
        kmer=ProfileMostProbableKmer(text, k, profile)
        Motifs.insert(i-1,kmer)
        if Score(Motifs)<Score(BestMotifs):
            BestMotifs=Motifs
    return BestMotifs

def main():
    """
    Demonstrates the usage of all major functions in this module with example data.
    """
    print("--- PatternCount ---")
    print(PatternCount("GCGCG", "GCG"))  # Output: 2

    print("\n--- FrequentWords ---")
    print(FrequentWords("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4))

    print("\n--- FrequencyMap ---")
    print(FrequencyMap("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4))

    print("\n--- PatternMatching ---")
    print(PatternMatching("ATAT", "GATATATGCATATACTT"))

    print("\n--- ReverseComplement ---")
    print(ReverseComplement("AAAACCCGGT"))

    print("\n--- SymbolArray & FasterSymbolArray ---")
    genome = "AAAAGGGGCCCCAAAAGGGGCCCC"
    print(SymbolArray(genome, "A"))
    print(FasterSymbolArray(genome, "A"))

    print("\n--- MinimumSkew & SkewArray ---")
    genome2 = "CATGGGCATCGGCCATACGCC"
    print(SkewArray(genome2))
    print(MinimumSkew(genome2))

    print("\n--- ApproximatePatternMatching & ApproximatePatternCount ---")
    print(ApproximatePatternMatching("CGCCCGAATCCAGAACGCATTCCCATATTTCATTGTAA", "ATTCTGGA", 3))
    print(ApproximatePatternCount("GAGG", "TTTAGAGCCTTCAGAGG", 2))

    print("\n--- HammingDistance ---")
    print(HammingDistance("GGGCCGTTGGT", "GGACCGTTGAC"))

    print("\n--- Count_old, Profile_old, Consensus, Score ---")
    motifs = ["ATCCAGCT", "GGGCAACT", "ATGGATCT", "AAGCAACC", "TTGGAACT", "ATGCCATT", "ATGGCACT"]
    print(Count_old(motifs))
    print(Profile_old(motifs))
    print(Consensus(motifs))
    print(Score(motifs))

    print("\n--- Pr, ProfileMostProbableKmer ---")
    profile = {'A': [0.2, 0.2, 0.3, 0.2, 0.3], 'C': [0.4, 0.3, 0.1, 0.5, 0.1], 'G': [0.3, 0.3, 0.5, 0.2, 0.4], 'T': [0.1, 0.2, 0.1, 0.1, 0.2]}
    print(Pr("CCGAG", profile))
    print(ProfileMostProbableKmer("ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT", 5, profile))

    print("\n--- GreedyMotifSearch ---")
    Dna = ["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"]
    print(GreedyMotifSearch(Dna, 3, 5))

    print("\n--- Count, Profile, RandomMotifs ---")
    print(Count(motifs))
    print(Profile(motifs))
    print(RandomMotifs(Dna, 3, 5))

    print("\n--- RandomizedMotifSearch ---")
    print(RandomizedMotifSearch(Dna, 3, 5))

    print("\n--- Motifs ---")
    pf = Profile(motifs)
    print(Motifs(pf, Dna))

    print("\n--- Normalize, WeightedDie, ProfileGeneratedString ---")
    probs = {'A': 0.1, 'C': 0.2, 'G': 0.3, 'T': 0.4}
    print(Normalize(probs))
    print(WeightedDie(probs))
    print(ProfileGeneratedString("ACGTACGT", pf, 3))

    print("\n--- GibbsSampler ---")
    print(GibbsSampler(Dna, 3, 5, 10))


if __name__ == "__main__":
    main()