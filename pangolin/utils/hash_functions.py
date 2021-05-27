#!/usr/bin/env python3
import hashlib
from Bio import SeqIO

def get_hash_string(record):
    seq = str(record.seq).upper().encode()
    hash_object = hashlib.md5(seq)
    hash_string = hash_object.hexdigest()
    return hash_string