#!/usr/bin/env python

import csv
from Bio import SeqIO
import os
import gzip
from pangolin.utils.log_colours import green,cyan
from pangolin.utils.report_collation import usher_parsing
from pangolin.utils.config import *

##### Report options #####


##### Target rules #####

rule all:
    input:
        os.path.join(config[KEY_TEMPDIR],"cache_assigned.csv"),
        os.path.join(config[KEY_TEMPDIR],"inference_report.csv")


rule usher_cache:
    input:
        fasta = os.path.join(config[KEY_TEMPDIR],"pass_qc.fasta"),
        hash_map = os.path.join(config[KEY_TEMPDIR],"hash_map.csv")
    output:
        cached = os.path.join(config[KEY_TEMPDIR],"cache_assigned.csv"),
        for_inference = os.path.join(config[KEY_TEMPDIR],"not_cache_assigned.fasta")
    params:
        cache = config[KEY_ASSIGNMENT_CACHE]
    run:
        if not params.cache:
            # Copy input.fasta to output.for_inference (all sequences need inference)
            with open(output.for_inference, "w") as fseq:
                for record in SeqIO.parse(input.fasta, "fasta"):
                    fseq.write(f">{record.description}\n{record.seq}\n")
            # Make a header-only CSV output, just like FINAL_HEADER but with hash instead of taxon:
            with open(output.cached,"w") as fw:
                cache_header = FINAL_HEADER.copy()
                cache_header[0] = "hash"
                fw.write(",".join(cache_header) + '\n')
        else:
            # Make a set of input sequence hashes; the hashes that are found in the cache will
            # be deleted so that after checking the cache, the leftovers will be the set of
            # sequence hashes (used as names in input.fasta) that need inference.
            seqs_to_assign = set()
            with open(input.hash_map, "r") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    seqs_to_assign.add(row["hash"])

            # Scan through cache to find hashes in seqs_to_assign.
            with gzip.open(params.cache,"rt") as f, \
                open(output.cached,"w") as fw:
                reader = csv.DictReader(f)
                fieldnames = reader.fieldnames[:]
                writer = csv.DictWriter(fw, fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
                writer.writeheader()
                for row in reader:
                    hash = row["hash"]
                    if hash in seqs_to_assign:
                        cache_note = "Assigned from cache"
                        if not row["note"]:
                            row["note"] = cache_note
                        elif cache_note not in row["note"]:
                            row["note"] += "; " + cache_note
                        writer.writerow(row)
                        # Remove this cached hash from seqs_to_assign.
                        seqs_to_assign.remove(hash)

            # Sequences left over in seqs_to_assign were not found in cache and will need inference.
            with open(output.for_inference, "w") as fseq:
                for record in SeqIO.parse(input.fasta, "fasta"):
                    if record.id in seqs_to_assign:
                        fseq.write(f">{record.id}\n{record.seq}\n")


rule usher_inference:
    input:
        fasta = rules.usher_cache.output.for_inference,
        reference = config[KEY_REFERENCE_FASTA],
        usher_protobuf = config[KEY_USHER_PB]
    params:
        vcf = os.path.join(config[KEY_TEMPDIR], "sequences.aln.vcf")
    threads: workflow.cores
    output:
        txt = os.path.join(config[KEY_TEMPDIR], "clades.txt")
    log:
        os.path.join(config[KEY_TEMPDIR], "logs/usher.log")
    shell:
        """
        echo "Using UShER as inference engine."
        if [ -s {input.fasta:q} ]; then
            faToVcf -includeNoAltN <(cat {input.reference:q} <(echo "") {input.fasta:q}) {params.vcf:q}
            usher -n -D -i {input.usher_protobuf:q} -v {params.vcf:q} -T {workflow.cores} -d '{config[tempdir]}' &> {log}
        else
            rm -f {output.txt:q}
            touch {output.txt:q}
        fi
        """

rule usher_to_report:
    input:
        txt = rules.usher_inference.output.txt,
    output:
        csv = os.path.join(config[KEY_TEMPDIR],"inference_report.csv")
    run:
        usher_parsing(input.txt, output.csv)
