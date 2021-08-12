import argparse
from python.common import main

def parse_input():
    ## Just take in the fasta files and then the output name
    ## use the default path for the data

    purpose = ''' This is a scipt to use the Pasteur MLST profiles and define these for a set of input 
    acinetbacter baumannii isolates. Usage:
    python acba_mlst_runner.py --seqs <list_of_fastas>  --output <output_prefixes> --threads <num_cores_to_use>'''

    parser = argparse.ArgumentParser(description=purpose,
                                     prog='acba_mlst_runner.py')

    parser.add_argument('--seqs', required=True, help='List of seqeuence files (FASTA) (required)', type=str)
    parser.add_argument('--output', required=True, help='Prefix of output files  (required)', type=str)
    parser.add_argument('--threads', default=1, help='Number of threads to use for ORF finder', type = int)

    args = parser.parse_args()

    return args

def main_run():
    args = parse_input()
    main(args)


if __name__ == '__main__':
    main_run()

