import re
import sys
import pandas
import time
import os

import pandas as pd

from python import orfipy_search
from python import blast_search
from python import blast_processing


def main(input_args):
    tic_setup = time.perf_counter()
    print("Beginning setup")
    pandas.set_option('display.max_columns', 500)

    # parse command line


    ###############################################################################
    ## Lets go through the gff list in a for loop, we'll find the position of the #
    ## gene of interest and then from there extract this sequence and compare it ##
    ## to the gene alignment to test which is the right gene ######################
    ###############################################################################
    seq_files = open(input_args.seqs, "r")
    seq_lines = seq_files.read().splitlines()

    python_dir_name = os.path.dirname(os.path.realpath(__file__))
    data_dir = re.sub("python","data/", python_dir_name)
    print(data_dir)
    aa_dir_name = "./" +  input_args.output + "_aa_dir" + "/"
    if not os.path.exists(aa_dir_name):
        os.mkdir(aa_dir_name)

    df_names = ['id','ST','cpn60','fusA','gltA','pyrG','recA','rplB','rpoB','clonal_complex','species']
    out_df = pd.DataFrame(index=range(len(seq_lines)), columns =df_names)

    skip = False
    toc_setup = time.perf_counter()
    print("Took this long for initial set up: %s" % (toc_setup - tic_setup))
    print("Beginning iso run")
    for k, fasta_file in enumerate(seq_lines):
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("On isolate %s of %s" % (k + 1, len(seq_lines)))
        tic_iso_run = time.perf_counter()
        bassio_nameo = os.path.basename(fasta_file)
        print(bassio_nameo)

        ## Perform the ORF finder for the isolate

        #ORF_loc = orfipy_search.orfipy_search(fasta_file, input_args.threads)

        ## Search for cpn60
        cpn60_seq = blast_search.blast_search_for_gene(fasta_file, "cpn60", data_dir, aa_dir_name, bassio_nameo, input_args.threads)

        # Search for fusA
        fusA_seq = blast_search.blast_search_for_gene(fasta_file, "fusA", data_dir, aa_dir_name, bassio_nameo, input_args.threads)

        # Search for gltA
        gltA_seq = blast_search.blast_search_for_gene(fasta_file, "gltA", data_dir, aa_dir_name, bassio_nameo, input_args.threads)
        # Search for pyrG
        pyrG_seq = blast_search.blast_search_for_gene(fasta_file, "pyrG", data_dir, aa_dir_name, bassio_nameo, input_args.threads)
        # Search for recA
        recA_seq = blast_search.blast_search_for_gene(fasta_file, "recA", data_dir, aa_dir_name, bassio_nameo, input_args.threads)
        # Search for rplB
        rplB_seq = blast_search.blast_search_for_gene(fasta_file, "rplB", data_dir, aa_dir_name, bassio_nameo, input_args.threads)
        # Search for rpoB
        rpoB_seq = blast_search.blast_search_for_gene(fasta_file, "rpoB", data_dir, aa_dir_name, bassio_nameo, input_args.threads)
        toc_blast_run = time.perf_counter()
        print("Took this long for isolate blast searching %s (s)" % (toc_blast_run - tic_iso_run))
        tic_blast_process = time.perf_counter()

        res_dict = { "id" : bassio_nameo,
                   "cpn60" : blast_processing.process_blast(cpn60_seq, 405,"cpn60"),
                     "fusA" : blast_processing.process_blast(fusA_seq, 633,"fusA"),
                    "gltA" : blast_processing.process_blast(gltA_seq, 483, "gltA"),
                     "pyrG" : blast_processing.process_blast(pyrG_seq, 297, "pyrG"),
                     "recA" : blast_processing.process_blast(recA_seq, 372, "recA"),
                     "rplB" : blast_processing.process_blast(rplB_seq, 330, "rplB"),
                     "rpoB" : blast_processing.process_blast(rpoB_seq, 456, "rpoB")}
        res_df = pd.Series(res_dict).to_frame().transpose()
        res_df.to_csv("./example_res_df.csv", index=None)

        out_st = blast_processing.ST_process(data_dir, res_df)
        out_df.loc[k] = out_st.iloc[0]
        toc_blast_run = time.perf_counter()
        print("Took this long for isolate blast processing %s (s)" % (toc_blast_run - tic_blast_process))

    toc_blast_run = time.perf_counter()
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print()
    print("Took this long for %s isolates %s (minutes)" % (len(seq_lines),round((toc_blast_run - tic_setup)/ 60, 3)))
    out_df.to_csv((input_args.output + "_ST.csv"), index=None)








