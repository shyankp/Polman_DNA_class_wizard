#!/usr/bin/env python
# coding: utf-8


# # Helpful Variables
standard_code = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "UGA": "*",
    "UGU": "C",
    "UGC": "C",
    "UGG": "W",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


aa_mol_weights = {
    "A": 89.09,
    "C": 121.15,
    "D": 133.1,
    "E": 147.13,
    "F": 165.19,
    "G": 75.07,
    "H": 155.16,
    "I": 131.17,
    "K": 146.19,
    "L": 131.17,
    "M": 149.21,
    "N": 132.12,
    "P": 115.13,
    "Q": 146.15,
    "R": 174.2,
    "S": 105.09,
    "T": 119.12,
    "V": 117.15,
    "W": 204.23,
    "X": 0,
    "Y": 181.19,
}


# # Class Seq


class seq:
    # Call instance attributes
    def __init__(self, name, type, organism, sequence):
        self.name = name
        self.type = type
        self.organism = organism
        self.sequence = sequence

    # define the info function
    def info(self):
        print(self.name)
        print(self.type)
        print(self.organism)
        print(self.sequence)

    # define the length function
    def length(self):
        x = len(self.sequence)
        print(x)

    # define the fasta_out function.
    # write to a file using the sequence name as part of the file name

    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            ">"
            + self.name
            + "_"
            + self.organism
            + "_"
            + self.type
            + "\n"
            + self.sequence
        )
        f.close()


# # Class Protein


class protein(seq):
    def __init__(self, name, type, organism, sequence, size):
        super().__init__(name, type, organism, sequence)
        self.size = size

    def info(self):
        print(self.name)
        print(self.type)
        print(self.organism)
        print(self.sequence)
        print(self.size)

    def mol_weight(self):
        total_weight = 0
        for aa in self.sequence:
            weight = aa_mol_weights.get(aa, 0)
            total_weight += weight
        return total_weight

    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            ">"
            + self.name
            + "_"
            + self.organism
            + "_"
            + self.type
            + "_"
            + self.size
            + "\n"
            + self.sequence
        )
        f.close()


# # Class Nucleotide


class nucleotide(seq):
    def __init__(self, name, type, organism, sequence):
        super().__init__(name, type, organism, sequence)

    def gc_content(self):
        gc_count = self.sequence.count("G") + self.sequence.count("C")
        total_base = len(self.sequence)
        gc_percentage = (gc_count / total_base) * 100
        print(gc_percentage)


# # Class DNA


class DNA(nucleotide):
    def __init__(self, name, type, organism, sequence):
        super().__init__(name, type, organism, sequence)

    # For the transcribe method
    def transcribe(self):
        rna_seq = self.sequence.replace("T", "U")
        print(rna_seq)

    # make the reverse completment
    def reverse_complement(self):
        complement_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}
        reverse_sequence = self.sequence[::-1]
        rev_seq = ""
        for base in reverse_sequence:
            complement_base = complement_dict.get(base, base)
            rev_seq += complement_base
        return rev_seq

    def six_frames(self):
        reverse = self.reverse_complement()
        frames = []
        # fwd frames
        for i in range(3):
            frames.append(self.sequence[i:])
        # rev frames
        for i in range(3):
            frames.append(reverse[i:])
        return frames


# # Class RNA


class RNA(nucleotide):
    def __init__(self, name, type, organism, sequence):
        super().__init__(name, type, organism, sequence)

    def start(self):
        start_codon = "AUG"
        index = self.sequence.find(start_codon)
        return index

    def translate(self):
        start_index = self.start()
        if start_index == -1:
            return "No start codon"

        cl = []
        for i in range(start_index, len(self.sequence), 3):
            codon = self.sequence[i : i + 3]
            if len(codon) == 3:
                cl.append(codon)

        aa_sequence = ""
        for codon in cl:
            aa_sequence += standard_code.get(codon)
        return aa_sequence


# # Test the new methods


uidA = DNA(
    name="uidA_DNA",
    organism="bacteria",
    type="DNA",
    sequence="CGCATGTTACGTCCTGTAGAAACCCCAACCCGTGAAATCAAAAAA",
)

uidA.fasta_out()


uidA6 = uidA.six_frames()
print(uidA6)

uidA_rev = uidA.reverse_complement()
print(uidA_rev)


uidA.transcribe()


uidA_RNA = RNA(
    name="uidA_RNA",
    organism="bacteria",
    type="RNA",
    sequence="CGCAUGUUACGUCCUGUAGAAACCCCAACCCGUGAAAUCAAAAAA",
)

uidA_RNA.fasta_out()


uidA_translate = uidA_RNA.translate()
print(uidA_translate)


uidA_protein = protein(
    name="uidA_protein",
    organism="bacteria",
    type="protein",
    sequence="MLRPVETPTREIKK",
    size="30",
)

uidA_protein.fasta_out()
uidA_mw = uidA_protein.mol_weight()
print(uidA_mw)
