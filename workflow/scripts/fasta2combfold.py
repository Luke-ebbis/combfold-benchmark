#!/usr/bin/env python
# -*- coding: utf-8 -*-
import Bio.PDB
from Bio.SeqIO import FastaIO
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
    tot = 0
    for record in records:
        current = dict()
        chainid = record.name.replace(":", "_").split("_")[1]
        name = f"{chainid}0"
        current['name'] = name
        current['chain_names'] = [chainid]
        current['start_res']= 1
        seq = str(record.seq)
        current['sequence'] = seq
        if len(seq) > 1800:
            raise RuntimeError(f"Sequence is {len(seq)}. This is too long")
        current['sequence'] = seq
        tot += len(seq)
        logging.info(f"Found record {name} with {chainid} -- len: {len(seq)}")
        if "X" not  in seq:
            job[name] = current
        else:
            current['sequence'] = "".join([s for s in seq if s != "X"])
            logger.info(f"In {name} X residues are removed")
            job[name] = current
        logger.info(f"Total amount of residues: {tot}")
    return job


def process_combfold_dict(records):
    """Remove duplicated sequences and put them under the same name 
    """
    updated_records = records.copy()
    removed_names = set()
    for name in records:
        logging.debug(f"checking {name} for duplicated sequence chains")
        current_sequence = records[name]['sequence']
        for check_name in records:
            if check_name != name:
                current_sequence = records[name]['sequence']
                logging.debug(f"checking {check_name} ")
                other_sequence = records[check_name]['sequence']
                if other_sequence == current_sequence and name not in removed_names:
                    removed = updated_records.pop(check_name)
                    removed_names.add(check_name)
                    removed_chain_names = [l for l in removed['chain_names']]
                    logging.debug(f"Sequence duplicate in {name}-{check_name}"\
                                  f"\nChain {removed_chain_names} will be added")
                    updated_records[name]['chain_names'] += removed_chain_names
    return updated_records


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
    processed_job = process_combfold_dict(job)
    write_job(processed_job, args.output)


if __name__ == "__main__":
    main()
