import os

rule all:
    input:
        os.path.join(config["tempdir"], "classified.csv")

rule minimap2_to_reference:
    input:
        fasta = config["query"],
        reference = config["reference"]
    output:
        sam = os.path.join(config["tempdir"], "query_mapped.sam")
    shell:
        """
        minimap2 -ax asm5 {input.reference} {input.fasta} > {output}
        """

rule find_query_snps:
    input:
        sam = rules.minimap2_to_reference.output.sam,
        ref = config["reference"]
    output:
        snps = os.path.join(config["tempdir"], "query_snps.csv")
    shell:
        """
        sam_2_snps.py -s {input.sam} -r {input.ref} -o {output.snps}
        """

rule classify:
    input:
        query_snps = rules.find_query_snps.output.snps,
        defining_snps = config["defining_snps"]
    output:
        csv = os.path.join(config["tempdir"], "classified.csv")
    shell:
        """
        snp_based_classifier.py -snps {input.defining_snps} -q {input.query_snps} -o {output.csv}
        """
