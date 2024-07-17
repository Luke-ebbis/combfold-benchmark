#!/usr/bin/env python
# -*- coding: utf-8 -*-
import Bio.PDB
from Bio.SeqIO import PdbIO, FastaIO
from Bio import SeqIO
import os
import itertools
import argparse
from typing import List, Tuple, Set, Dict
import logging
import json
import sys
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    filemode='w')
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.DEBUG)
console_handler.setFormatter(
logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger = logging.getLogger()


def read_fasta(input_file):
    records = list(SeqIO.parse(input_file, "fasta"))
    logging.info(f"Read in {len(records)} chains from {input_file}")
    return records


def format_combfold_dict(records):
    # TODO Handle X amino acid. Handle concat of unique seqs
    job = dict()
    for record in records:
        current = dict()
        chainid = record.name.replace(":", "_").split("_")[1]
        name = f"{chainid}0"
        current['name'] = name
        current['chain_names'] = [chainid]
        current['start_res']= 1
        seq = str(record.seq)
        current['sequence'] = seq
        job[name] = current
        if "X" in seq:
            continue
        if len(seq) > 1800:
            raise RuntimeError(f"Sequence is {len(seq)}. This is too long")
        current['sequence'] = seq
        logging.info(f"Found record {name} with {chainid}")
        job[name] = current
    return job


def write_job(jobdict: Dict, file) -> None:
    with open(file, "w") as f:
        json.dump(jobdict, f)
        logging.info(f"Wrote job to {file}")

def main():
    parser = argparse.ArgumentParser(
                        prog='Combfoldjobinator',
                        description=f'Extract the primary structure from a protein file',
                        epilog='Sibbe Bakker')
    parser.add_argument('fasta',
                        help="fasta file.")
    parser.add_argument('output',
                        help="Name of the combfold job config")
    args = parser.parse_args()

    sequences = read_fasta(args.fasta)
    job = format_combfold_dict(sequences)
    write_job(job, args.output)


if __name__ == "__main__":
    main()
