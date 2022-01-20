import gzip
import os
import sys



def sequence_qc(in_fasta,out_qc_pass,out_qc_fail):
    print(green("****\nRunning sequence QC"))

    if os.path.exists(in_fasta):
        file_ending = query.split(".")[-1]
        if file_ending in ["gz","gzip","tgz"]:
            query = gzip.open(query, 'rt')
        elif file_ending in ["xz","lzma"]:
            query = lzma.open(query, 'rt')
            
    post_qc_query = os.path.join(tempdir, 'query.post_qc.fasta')
    fw_pass = open(post_qc_query,"w")
    qc_fail = os.path.join(tempdir,'query.failed_qc.fasta')
    fw_fail = open(qc_fail,"w")

    total_input = 0
    total_pass = 0
    
    try:
        for record in SeqIO.parse(query, "fasta"):
            total_input +=1
            record.description = record.description.replace(' ', '_').replace(",","_")
            record.id = record.description
            if "," in record.id:
                record.id=record.id.replace(",","_")

            if len(record) <args.minlen:
                record.description = record.description + f" fail=seq_len:{len(record)}"
                fw_fail.write(f">{record.description}\n{record.seq}\n")
            else:
                num_N = str(record.seq).upper().count("N")
                prop_N = round((num_N)/len(record.seq), 2)
                if prop_N > args.maxambig:
                    record.description = record.description + f" fail=N_content:{prop_N}"
                    fw_fail.write(f">{record.description}\n{record.seq}\n")
                else:
                    total_pass +=1
                    seq = str(record.seq).replace("-","")
                    fw_pass.write(f">{record.description}\n{seq}\n")
    except UnicodeDecodeError:
        sys.stderr.write(cyan(
            f'Error: the input query fasta could not be parsed.\n' +
            'Double check your query fasta and that compressed stdin was not passed.\n' +
            'Please enter your fasta sequence file and refer to pangolin usage at: https://cov-lineages.org/pangolin.html' +
            ' for detailed instructions.\n'))
        sys.exit(-1)

    print(green("Number of sequences detected: ") + f"{total_input}")
    print(green("Total passing QC: ") + f"{total_pass}")
    fw_fail.close()
    fw_pass.close()

    if total_pass == 0:
        with open(outfile, "w") as fw:
            fw.write("taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,version,pangolin_version,pangoLEARN_version,pango_version,status,note\n")
            for record in SeqIO.parse(os.path.join(tempdir,'query.failed_qc.fasta'), "fasta"):
                desc = record.description.split(" ")
                reason = ""
                for item in desc:
                    if item.startswith("fail="):
                        reason = item.split("=")[1]
                fw.write(f"{record.id},None,,,,,,PANGO-{PANGO_VERSION},{__version__},{pangoLEARN.__version__},{PANGO_VERSION},fail,{reason}\n")
        print(cyan(f'Note: no query sequences have passed the qc\n'))
        sys.exit(0)