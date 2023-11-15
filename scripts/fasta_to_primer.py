import json
from Bio import SeqUtils
from Bio.Seq import Seq
from Bio import SeqIO  
import primer3

def load_design_parameters(config_file):
    with open(config_file, 'r') as f:
        design_params = json.load(f)
    return design_params

def calculate_tm(primer_sequence):
    result = primer3.calc_tm(primer_sequence)
    return result

def design_primers(sequence, design_params):
    seq = Seq(sequence)
    primers = []
    #It should only be 5 primers
    count = 5 

    min_length = design_params['min_length']
    max_length = design_params['max_length']
    min_gc = design_params['min_gc'] / 100
    max_gc = design_params['max_gc'] / 100
    max_tm = design_params['max_tm']

    for i in range(len(seq) - min_length + 1):
        for j in range(i + min_length, min(i + max_length, len(seq)) + 1):
            primer = seq[i:j]

            gc_content = SeqUtils.gc_fraction(primer)

            tm = calculate_tm(str(primer))

            # Check if criteria are met
            if min_gc <= gc_content <= max_gc and tm <= max_tm:
                primers.append(primer)
                count -= 1
                if count == 0:
                    return primers

    return primers

if __name__ == '__main__':
    fasta_file = 'gene.fasta'
    config_file = 'primer_design_config.json'

    with open(fasta_file, 'r') as seq_file:
        seq_record = SeqIO.read(seq_file, "fasta")
        sequence = str(seq_record.seq)

    design_params = load_design_parameters(config_file)

    designed_primers = design_primers(sequence, design_params)

    # Print or save the designed primers
    print("Designed Primers:")
    for primer in designed_primers:
        print(primer)
