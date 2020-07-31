# DNA Toolkit file
from structures import *
from collections import Counter


def validate_seq(seq):
    """Check the sequence to make sure it is a valid DNA string"""

    tmpseq = seq.upper()
    for nuc in tmpseq:
        if nuc not in DNA_Nucleotides:
            return False
    return True


def nucleotide_frequency(seq):
    """Count nucleotides in a given sequence. Return a dictionary"""

    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict

    # More Pythonic, using Counter
    # return dict(Counter(seq))


def transcription(seq):
    """DNA -> RNA Transcription. Replacing Thymine with Uracil"""
    return seq.replace("T", "U")


def reverse_complement(seq):
    """
    Swapping adenine with thymine and guanine with cytosine.
    Reversing newly generated string
    """
    return ''.join([DNA_ReverseComplement[nuc] for nuc in seq])[::-1]

    # Pythonic approach. A little bit faster solution.
    # mapping = str.maketrans('ATCG', 'TAGC')
    # return seq.translate(mapping)[::-1]


def gc_content(seq):
    """GC Content in a DNA/RNA sequence"""
    return round((seq.count('C') + seq.count('G')) / len(seq) * 100)


def gc_content_subsec(seq, k=20):
    """GC Content in a DNA/RNA sub-sequence length k. k=20 by default"""
    res = []
    for i in range(0, len(seq) - k + 1, k):
        subseq = seq[i:i + k]
        res.append(gc_content(subseq))
    return res


def translate_seq(seq, init_pos=0):
    """Translates a DNA sequence into an aminoacid sequence"""
    return [DNA_Codons[seq[pos:pos + 3]] for pos in range(init_pos, len(seq) - 2, 3)]


def codon_usage(seq, aminoacid):
    """Provides the frequency of each codon encoding a given aminoacid in a DNA sequence"""
    tmpList = []
    for i in range(0, len(seq) - 2, 3):
        if DNA_Codons[seq[i:i + 3]] == aminoacid:
            tmpList.append(seq[i:i + 3])

    freqDict = dict(Counter(tmpList))
    totalWight = sum(freqDict.values())
    for seq in freqDict:
        freqDict[seq] = round(freqDict[seq] / totalWight, 2)
    return freqDict

def gen_reading_frames(seq):
    """Generate the six reading frames of a DNA sequence, including reverse complement"""
    frames = []
    for i in range(3):
        frames.append(translate_seq(seq, i))
        frames.append(translate_seq(reverse_complement(seq), i))
    return frames

def proteins_from_rf(aa_seq):
    """Compute all possible proteins in an aminoacid seq and return a list of possible proteins"""
    current_prot = []
    proteins = []
    queue = []
    for aa in aa_seq:
        if len(queue) >= 6:
            queue.pop(0)
        queue.append(aa)
        if aa == "*":
            # STOP accumulating amino acids if _ - STOP was found
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            # START accumulating amino acids if M - START was found
            if aa == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa

    return proteins


#def slippage_detect(seven_amino):
#    """it must be N NNW WWZ
#    (codons are shown in the incoming or 0-frame),
#    where N is a stretch of three identical
#    nucleotides, W is either AAA or UUU, and Z ≠
#    G."""
#
#    (AAA|TTT|GGG|CCC)(AAA|UUU)(A|T|C)