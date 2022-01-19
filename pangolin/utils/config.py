
KEY_QUERY_FASTA="query_fasta"
KEY_OUTDIR="outdir"
KEY_OUTFILE="outfile"
KEY_TEMPDIR="tempdir"
KEY_ALIGNDIR="aligndir"
KEY_ALIGNMENT_OUT="alignment_out"

KEY_ANALYSIS_MODE="analysis_mode"

KEY_TRIM_START="trim_start"
KEY_TRIM_END="trim_end"
KEY_QC_FAIL="qc_fail"
KEY_ALIAS_FILE="alias_file"
KEY_CONSTELLATION_FILES="constellation_files"
KEY_SKIP_DESIGNATION_HASH="skip_designation_hash"
KEY_USE_CACHE="use_cache"
KEY_VERBOSE="verbose"

KEY_PANGOLEARN_VERSION="pangoLEARN_version"
KEY_PANGOLIN_VERSION="pangolin_version"
KEY_CONSTELLATIONS_VERSION="constellation_version"
KEY_SCORPIO_VERSION="scorpio_version"
KEY_PANGO_VERSION="pango_version"

KEY_THREADS="threads"

dependency_list = ["gofasta","minimap2","snakemake"]
module_list = ["Bio","sklearn","pandas","joblib","pysam","pangoLEARN","constellations"]
