#!/usr/bin/env python3

from Bio import SeqIO
import pysam
import re, itertools, operator
import sys
import argparse

lambda_dict = {'M': (lambda query_start, ref_start, length, seq: (query_start + length, ref_start + length, seq[query_start:query_start + length] )),
               'I': (lambda query_start, ref_start, length, seq: (query_start + length, ref_start         , ''                                    )),
               'D': (lambda query_start, ref_start, length, seq: (query_start         , ref_start + length, '-' * length                          )),
               'N': (lambda query_start, ref_start, length, seq: (query_start         , ref_start + length, '-' * length                          )),
               'S': (lambda query_start, ref_start, length, seq: (query_start + length, ref_start         , ''                                    )),
               'H': (lambda query_start, ref_start, length, seq: (query_start         , ref_start         , ''                                    )),
               'P': (lambda query_start, ref_start, length, seq: (query_start         , ref_start         , ''                                    )),
               '=': (lambda query_start, ref_start, length, seq: (query_start + length, ref_start + length, seq[query_start:query_start + length] )),
               'X': (lambda query_start, ref_start, length, seq: (query_start + length, ref_start + length, seq[query_start:query_start + length] ))}

def parse_sam_line(AlignedSegment):
    """
    d is a dictionary with SAM field names as keys and their value from
    one line of a SAM alignment as values
    """
    line = str(AlignedSegment)
    names = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL']
    d = {x: y for x,y in zip(names, line.split()[0:11])}
    return(d)


def split_sam_cigar_operation(one_operation):
    """
    o is a tuple with the format(operation, size)
    e.g. ('M', 2377) or ('D', 1)
    """
    type = one_operation[-1:]
    size = int(one_operation[:-1])
    o = (type, size)
    return(o)


def split_sam_cigar(cigar):
    """
    m is a list of strings (e.g. ['1000M', '4I'])
    """
    r = re.compile('\d{1,}[A-Z]{1}')
    l = re.findall(r, cigar)
    return(l)


def get_sam_cigar_operations(cigar):
    """
    operations is a list of tuples that correspond to
    operations to apply in order:
    [('M', 10000),('I', 3),('M', 19763)]
    """
    operations_raw = split_sam_cigar(cigar)
    operations = [split_sam_cigar_operation(x) for x in operations_raw]
    return(operations)


def get_one_string(sam_line, rlen):
    """
    Transform one line of the SAM alignment into sample sequence in unpadded
    reference coordinates (insertions relative to the reference are omitted
    """

    # parsed sam line
    aln_info_dict = parse_sam_line(sam_line)

    # query sequence name
    QNAME = aln_info_dict['QNAME']

    # CIGAR STRING
    CIGAR = aln_info_dict['CIGAR']

    # According to the SAM spec:
    # "POS: 1-based leftmost mapping POSition of the first CIGAR operation that
    # “consumes” a reference base (see table above)."
    # But note that pysam converts this field to 0-bsaed coordinates for us
    POS = int(aln_info_dict['POS'])

    # Query seq:
    SEQ = aln_info_dict['SEQ']

    # According to the SAM spec:
    # "If POS < 1, unmapped read, no assumptions can be made about RNAME and CIGAR"
    # But note that pysam converts POS to 0-bsaed coordinates for us (as above)
    if POS < 0:
        # TO DO: SOME SENSIBLE RETURN HERE
        return(None)

    # parse the CIGAR string to get the operations:
    operations = get_sam_cigar_operations(CIGAR)

    # left-pad the new sequence with gaps if required
    new_seq = '*' * POS

    # then build the sequence:
    qstart = 0
    rstart = POS
    for o in operations:
        operation = o[0]
        size = o[1]

        # based on this CIGAR operation, call the relavent lambda function
        # from the dict of lambda functions, returns sequence to be appended
        # and the next set of coordinates
        new_qstart, new_rstart, extension = lambda_dict[operation](qstart, rstart, size, SEQ)

        new_seq = new_seq + extension

        qstart = new_qstart
        rstart = new_rstart

    rightpad = '*' * (rlen - len(new_seq))

    new_seq = new_seq + rightpad

    return(new_seq)


def check_and_get_flattened_site(site):
    """
    A per-site check that there isn't any ambiguity between
    alignments within a single sequence
    """

    check = sum([x.isalpha() for x in site])
    if check > 1:
        sys.stderr.write('ambiguous overlapping alignment')
        # sys.exit()
        return('&')

    # because {A, C, G, T} > {-} > {*}, we can use max()
    base = max(site)
    return(base)


def get_seq_from_block(sam_block, rlen):

    block_lines_sites_list = [get_one_string(sam_line, rlen) for sam_line in sam_block]

    if len(block_lines_sites_list) == 1:
        return(block_lines_sites_list[0])

    elif len(block_lines_sites_list) > 1:
        # # as an alternative to check_and_get_flattened_site() we can flatten
        # # the site with no checks (about three times as fast):
        # flattened_site_list = [max(x) for x in zip(*[list(x) for x in block_lines_sites_list])]

        flattened_site_list = [check_and_get_flattened_site(x) for x in zip(*[list(x) for x in block_lines_sites_list])]
        seq_flat = ''.join(flattened_site_list)

        # replace central '*'s with 'N's, and external '*'s with '-'s
        return(seq_flat)


def get_snp_list(reference, query):
    lr = list(reference)
    lq = list(query)

    snp_list = []
    for position, nucleotide_pair in enumerate(zip(lr, lq)):
        ref_nuc = nucleotide_pair[0].lower()
        que_nuc = nucleotide_pair[1].lower()

        if ref_nuc in ['a', 'c', 'g', 't'] and que_nuc in ['a', 'c', 'g', 't']:
            if ref_nuc != que_nuc:
                snp_list.append(str(position + 1) + ref_nuc.upper() + que_nuc.upper())

    return(snp_list)


def sam_2_snps(samfile, reference, output):

    samfile = pysam.AlignmentFile(samfile, 'r')
    reference = SeqIO.read(reference, 'fasta')

    # The length of the reference sequence:
    RLEN = samfile.header['SQ'][0]['LN']
    # The name of the target sequence:
    REF = samfile.header['SQ'][0]['SN']

    if REF != reference.id:
        sys.exit('reference names differ!')
    if RLEN != len(reference.seq):
        sys.exit('reference lengths differ!')

    with open(output, 'w') as file_out:
        for query_seq_name, one_querys_alignment_lines in itertools.groupby(samfile, lambda x: parse_sam_line(x)['QNAME']):
            # one_querys_alignment_lines is an iterator corresponding to all the lines
            # in the SAM file for one query sequence

            one_querys_alignment_lines = [x for x in one_querys_alignment_lines if parse_sam_line(x)['SEQ'] != 'None']
            if len(one_querys_alignment_lines) == 0:
                sys.stderr.write(query_seq_name + ' has 0-length SEQ field in alignment\n')
                continue

            seq = get_seq_from_block(sam_block = one_querys_alignment_lines, rlen = RLEN)

            snp_list = get_snp_list(reference.seq, seq)

            file_out.write(query_seq_name + ',' + ';'.join(snp_list) + '\n')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="turn a sam file into snp information")
    parser.add_argument("-r","--reference", help = "reference in fasta format")
    parser.add_argument("-s","--sam", help = "alignment in sam format")
    parser.add_argument("-o","--outfile", help = "csv file to write")
    args = parser.parse_args()

    sam_2_snps(samfile = args.sam, reference = args.reference, output = args.outfile)














#
