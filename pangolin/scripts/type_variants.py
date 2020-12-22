#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import sys

if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

def get_nuc_position_from_aa_description(cds, aa_pos):
    """
    given a CDS (eg. S) and the number of an amino acid in it, get the
    1-based start position of that codon in Wuhan-Hu-1 ref coordinates
    nuc_pos is an integer which is 1-based start pos of codon
    """

    # these coordinates are 1-based
    CDS_dict = {"orf1ab": ((266, 13468), (13468, 21555)),
                "orf1a":   (266, 13468),
                "orf1b":   (13468, 21555),
                "s":       (21563, 25384),
                "orf3a":   (25393, 26220),
                "e":       (26245, 26472),
                "m":       (26523, 27191),
                "orf6":    (27202, 27387),
                "orf7a":   (27394, 27759),
                "orf8":    (27894, 28259),
                "n":       (28274, 29533),
                "orf10" :  (29558, 29674)}

    if cds.lower() not in CDS_dict.keys():
        sys.stderr.write("I don't know about cds: %s \n" % cds)
        sys.stderr.write("please use one of: %s" % ",".join(CDS_dict.keys()))
        sys.exit(1)

    if cds.lower() == "orf1ab":
        if aa_pos <= 4401:
            parsed_cds = "orf1a"
        else:
            parsed_cds = "orf1b"
            aa_pos = aa_pos - 4401
    else:
        parsed_cds = cds

    cds_tuple = CDS_dict[parsed_cds.lower()]

    nuc_pos = cds_tuple[0] + ((aa_pos - 1) * 3)

    if nuc_pos > cds_tuple[1]:
        sys.stderr.write("invalid amino acid position for cds %s : %d" % (cds, aa_pos))
        sys.exit(1)

    return nuc_pos


def parse_variants_in(csvfilehandle, refseq):
    """
    read in a variants file and parse its contents and
    return something sensible.
    format of mutations csv is:
    snp:T6954C
    del:11288:9
    aa:orf1ab:T1001I
    returns variant_list which is a list of dicts of snps, aas and dels,
    one dict per variant. format of subdict varies by variant type
    """

    variant_list = []

    with open(csvfilehandle, "r") as f:
        for line in f:
            l = line.strip()
            lsplit = l.split(":")

            if lsplit[0] == "snp":
                type = lsplit[0]
                ref_allele = lsplit[1][0]
                ref_start = int(lsplit[1][1:-1])
                alt_allele = lsplit[1][-1]
                ref_allele_check = refseq[ref_start - 1]

                if ref_allele != ref_allele_check:
                    sys.stderr.write("variants file says reference nucleotide at position %d is %s, but reference sequence has %s" %(ref_start, ref_allele, ref_allele_check))
                    sys.exit(1)

                newsnprecord = {"type": type, "ref_start": ref_start, "ref_allele": ref_allele, "alt_allele": alt_allele}
                variant_list.append(newsnprecord)

            elif lsplit[0] == "aa":
                type = lsplit[0]
                cds = lsplit[1]
                ref_allele = lsplit[2][0]
                AA_pos = int(lsplit[2][1:-1])
                alt_allele = lsplit[2][-1]

                ref_start = get_nuc_position_from_aa_description(cds, AA_pos)
                ref_allele_check = refseq[ref_start - 1:ref_start + 2].translate()

                if ref_allele != ref_allele_check:
                    sys.stderr.write("variants file says reference amino acid in CDS %s at position %d is %s, but reference sequence has %s" %(cds, AA_pos, ref_allele, ref_allele_check))
                    sys.exit(1)

                newaarecord = {"type": type, "cds": cds, "ref_start": ref_start, "ref_allele": ref_allele, "alt_allele": alt_allele}
                variant_list.append(newaarecord)

            elif lsplit[0] == "del":
                length = int(lsplit[2])
                newdelrecord = {"type": lsplit[0], "ref_start": int(lsplit[1]), "length": length, "ref_allele": refseq[int(lsplit[1]) - 1:int(lsplit[1]) + length - 1]}
                variant_list.append(newdelrecord)

            else:
                sys.stderr.write("couldn't parse the following line in the config file: %s" % line)
                sys.exit()

    return(variant_list)


def type_variants(fasta_in, reference, variants_in, variants_out_handle, write_all_variants = False):
    reference_seq = list(SeqIO.parse(reference, "fasta"))[0].seq
    variant_list = parse_variants_in(variants_in, reference_seq)

    variants_out = open(variants_out_handle, "w")
    variants_out.write("query,ref_count,alt_count,other_count,fraction_alt")
    if write_all_variants:
        with open(variants_in, "r") as f:
            for line in f:
                variants_out.write(",")
                variants_out.write(line.strip())
    variants_out.write("\n")

    with open(fasta_in, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            ref_count = 0
            alt_count = 0
            oth_count = 0

            alleles_list = []

            for var in variant_list:
                if var["type"] == "snp":
                    query_allele = record.seq[var["ref_start"] - 1]
                    if query_allele == var["ref_allele"]:
                        ref_count += 1
                    elif query_allele == var["alt_allele"]:
                        alt_count += 1
                    else:
                        oth_count += 1

                    alleles_list.append(query_allele)


                if var["type"] == "aa":
                    try:
                        query_allele = record.seq[var["ref_start"] - 1:var["ref_start"] + 2].translate()
                    except:
                        oth_count += 1
                        alleles_list.append("X")
                        continue

                    if query_allele == var["ref_allele"]:
                        ref_count += 1
                    elif query_allele == var["alt_allele"]:
                        alt_count += 1
                    else:
                        oth_count += 1

                    alleles_list.append(str(query_allele))


                if var["type"] == "del":
                    query_allele = record.seq[var["ref_start"] - 1:var["ref_start"] + var["length"] - 1]
                    if query_allele == var["ref_allele"]:
                        ref_count += 1
                        alleles_list.append("ref")
                    elif query_allele == "-" * var["length"]:
                        alt_count += 1
                        alleles_list.append("del")
                    else:
                        oth_count += 1
                        alleles_list.append("X")


            variants_out.write(record.id + ",")
            variants_out.write(str(ref_count) + ",")
            variants_out.write(str(alt_count) + ",")
            variants_out.write(str(oth_count) + ",")
            variants_out.write(str(round(alt_count / (alt_count + ref_count + oth_count), 4)))

            if write_all_variants:
                variants_out.write(",")
                variants_out.write(",".join(alleles_list))

            variants_out.write("\n")

    variants_out.close()

    pass


def parse_args():
    parser = argparse.ArgumentParser(description="""type an alignment in Wuhan-Hu-1 coordinates for variants defined in a config file""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--fasta-in', dest = 'fasta_in', help='alignment to type, in fasta format')
    parser.add_argument('--variants-config', dest = 'variants_in', help="""config file containing variants to type. Format is like:
                                                    snp:T6954C
                                                    del:11288:9
                                                    aa:orf1ab:T1001I
                                                    snp and del positions are 1-based nucleotide position relative to the alignment
                                                    aa position is 1-based codon position relative to the cds
    No header line or comment lines are permitted""")
    parser.add_argument('--reference', help='Wuhan-Hu-1 in fasta format (for typing the reference allele at deletions and sanity checking the config file)')
    parser.add_argument('--variants-out', dest = 'variants_out', help='csv file to write')
    parser.add_argument('--append-genotypes', dest = 'append_genotypes', action = 'store_true', help='if invoked, write the genotype for each variant in the config file to the output')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = parse_args()
    type_variants(fasta_in = args.fasta_in,
                  reference = args.reference,
                  variants_in = args.variants_in,
                  variants_out_handle = args.variants_out,
                  write_all_variants = args.append_genotypes)
                  