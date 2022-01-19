KEY_ANALYSIS_MODE="analysis_mode"

KEY_SKIP_DESIGNATION_HASH="skip_designation_hash"
KEY_USE_CACHE="use_cache"

KEY_QUERY_FASTA="query_fasta"

KEY_OUTDIR="outdir"
KEY_OUTFILE="outfile"

KEY_ALIGNDIR="aligndir"
KEY_ALIGNMENT_FILE="alignment_file"
KEY_ALIGNMENT_OUT="alignment_out"

KEY_TEMPDIR="tempdir"
KEY_NO_TEMP = "no_temp"

KEY_DATADIR="datadir"

KEY_TRIM_START="trim_start"
KEY_TRIM_END="trim_end"

KEY_ALIAS_FILE="alias_file"

KEY_CONSTELLATION_FILES="constellation_files"

KEY_PANGOLEARN_VERSION="pangoLEARN_version"
KEY_PANGOLIN_VERSION="pangolin_version"
KEY_CONSTELLATIONS_VERSION="constellation_version"
KEY_SCORPIO_VERSION="scorpio_version"
KEY_PANGO_VERSION="pango_version"
KEY_PANGO_DESIGNATION_VERSION="pango_designation_version"

KEY_VERBOSE="verbose"
KEY_THREADS="threads"

dependency_list = ["gofasta","minimap2","snakemake"]
module_list = ["Bio","sklearn","pandas","joblib","pysam","pangoLEARN","constellations"]
