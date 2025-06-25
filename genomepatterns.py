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
    return Pattern[::-1]

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

def main():
    """
    Example usage of genome pattern functions.
    """
    genome = "ATGCAATCGACTACGATCAATCGAGGGCC" # Example genome string
    pattern = "ATC"                          # Pattern to search for
    k = 3                                    # Length of k-mers to consider

    print(f"Genome: {genome}")
    print(f"Pattern: {pattern}")
    print(f"k: {k}\n")

    # PatternCount
    count = PatternCount(genome, pattern)
    print(f"PatternCount: '{pattern}' appears {count} times in the genome.")

    # FrequencyMap
    freq_map = FrequencyMap(genome, k)
    print(f"FrequencyMap: {freq_map}")

    # FrequentWords
    frequent = FrequentWords(genome, k)
    print(f"FrequentWords: Most frequent {k}-mers are {frequent}")

    # PatternMatching
    positions = PatternMatching(pattern, genome)
    print(f"PatternMatching: '{pattern}' found at positions {positions}")

    # ReverseComplement
    rev_comp = ReverseComplement(pattern)
    print(f"ReverseComplement: Reverse complement of '{pattern}' is '{rev_comp}'")

    # Reverse
    rev = Reverse(pattern)
    print(f"Reverse: Reverse of '{pattern}' is '{rev}'")

    # Complement
    comp = Complement(pattern)
    print(f"Complement: Complement of '{pattern}' is '{comp}'")

if __name__ == "__main__":
    main()
