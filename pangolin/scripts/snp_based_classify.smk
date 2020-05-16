import os

rule minimap2_to_reference:
    input:
        fasta = config["query"],
        reference = config["reference"]
    output:
        sam = os.path.join(config["tempdir"], "query_mapped.sam")
    shell:
        """
        minimap2 -ax asm5 {input.reference} {input.fasta} -o {output}
        """

rule find_query_snps:
    input:
        sam = rules.minimap2_to_reference.output.sam,
        ref = config["reference"]
    output:
        snps = os.path.join(config["tempdir"], "query_snps.csv")
    shell:
        """
        python sam_2_snps.py -s {input.sam} -r {input.ref} -o {output.snps}

        ben script here
        take in .sam file and output a csv of
        query_id1,SNP1;SNP2;SNP3
        query_id2,SNP1;SNP2;SNP3

        to make the SNP1;SNP2;SNP3 bit I have been using this function:
        def snp_list_to_snp_string(snp_list):
            #turn a snp list into a `;`-separated string of snps that are sorted by
            #position in the genome
            snp_string = ";".join(sorted(snp_list, key = lambda x : int(x[:-2])))
            return snp_string
        """

rule classify:
    input:
        query_snps = rules.find_query_snps.output.snps,
        defining_snps = config["defining_snps"]
    output:
        os.path.join(config["tempdir"], "classified.csv")
    shell:
        """
        python snp_based_classifier.py -snps {input.defining_snps} -q {input.query_snps}
        """
