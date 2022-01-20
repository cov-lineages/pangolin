

rule cache_sequence_assign:
    input:
        fasta = rules.hash_sequence_assign.output.for_inference
    output:
        cached = os.path.join(config["tempdir"],"cache_assigned.csv"),
        for_inference = os.path.join(config["tempdir"],"not_cache_assigned.fasta")
    params:
        use_cache = config["use_cache"],
        cache = config["cache"]
    run:
        if not params.use_cache or params.cache == "":
            with open(output.for_inference, "w") as fseq:
                for record in SeqIO.parse(input.fasta, "fasta"):
                    fseq.write(f">{record.description}\n{record.seq}\n")
            with open(output.cached,"w") as fw:
                fw.write("taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,version,pangolin_version,pangoLEARN_version,pango_version,status,note\n")
        else:
            seq_hashes = {}
            token = ""
            for record in SeqIO.parse(input.fasta, "fasta"):
                if record.id!="reference":
                    hash_string = get_hash_string(record)
                    if hash_string in seq_hashes:
                        seq_hashes[hash_string].append(record.id)
                    else:
                        seq_hashes[hash_string] = [record.id]

            with gzip.open(params.cache,"rt") as f, \
                open(output.cached,"w") as fw:
                reader = csv.DictReader(f)
                fieldnames = reader.fieldnames[:]
                fieldnames.remove("hash")
                fieldnames.append("taxon")
                writer = csv.DictWriter(fw, fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
                writer.writeheader()
                for row in reader:
                    if row["hash"] in seq_hashes:
                        for seq_id in seq_hashes[row["hash"]]:
                            new_row = row
                            new_row["taxon"] = seq_id
                            del new_row["hash"]
                            if new_row["note"] == "":
                                new_row["note"] = "Assigned from cache"
                            else:
                                new_row["note"] += ";Assigned from cache"
                            writer.writerow(row)

            seqs_to_assign = set()
            for hash in seq_hashes:
                seqs_to_assign.update(seq_hashes[hash])

            with open(output.for_inference, "w") as fseq:
                for record in SeqIO.parse(input.fasta, "fasta"):
                    if record.id in seqs_to_assign:
                        fseq.write(f">{record.description}\n{record.seq}\n")
