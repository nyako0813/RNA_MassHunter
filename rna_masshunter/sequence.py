def normalize_rna_sequence(sequence: str) -> str:
    return sequence.upper().replace("T", "U").replace(" ", "").replace("\n", "")
