import sys
import json
from Bio import SeqUtils
from Bio.Seq import Seq
import primer3

def load_design_parameters(config_file):
    with open(config_file, 'r') as f:
        design_params = json.load(f)
    return design_params

def calculate_tm(primer_sequence):
    # Use primer3 to calculate Tm
    result = primer3.calcTm(primer_sequence)
    return result.tm

def design_primers(sequence, design_params):
    seq = Seq(sequence)
    primers = []

    min_length = design_params.get('min_length', 18)
    max_length = design_params.get('max_length', 24)
    min_gc = design_params.get('min_gc', 50)
    max_gc = design_params.get('max_gc', 60)
    max_tm = design_params.get('max_tm', 67)

    for i in range(len(seq) - min_length + 1):
        for j in range(i + min_length, min(i + max_length, len(seq)) + 1):
            primer = seq[i:j]

            # Calculate GC content
            gc_content = SeqUtils.GC(primer)

            # Calculate melting temperature (Tm)
            tm = calculate_tm(str(primer))

            # Check if criteria are met
            if min_gc <= gc_content <= max_gc and tm <= max_tm:
                primers.append(primer)

    return primers

def calculate_tm(primer_sequence):
    # Use primer3 to calculate Tm
    result = primer3.calcTm_Wallace(primer_sequence)
    return result.tm


if __name__ == '__main__':
    config_file = 'primer_design_config.json'  
    sequence = 'gene.fasta' 

    design_params = load_design_parameters(config_file)
    designed_primers = design_primers(sequence, design_params)

    # Print or save the designed primers
    for primer in designed_primers:
        print(primer)
