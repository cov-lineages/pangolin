KEY_ANALYSIS_MODE="analysis_mode"

KEY_QUERY_FASTA="query_fasta"
KEY_REFERENCE_FASTA="reference_fasta"

KEY_OUTDIR="outdir"
KEY_OUTFILE="outfile"

KEY_ALIGNDIR="aligndir"
KEY_ALIGNMENT_FILE="alignment_file"
KEY_ALIGNMENT_OUT="alignment_out"

KEY_TEMPDIR="tempdir"
KEY_NO_TEMP = "no_temp"

KEY_DATADIR="datadir"

KEY_MINLEN="minlen"
KEY_MAXAMBIG="maxambig"
KEY_TRIM_START="trim_start"
KEY_TRIM_END="trim_end"

KEY_ALIAS_FILE="alias_file"

KEY_CONSTELLATION_FILES="constellation_files"
KEY_USHER_PB = "usher_pb"
KEY_PLEARN_MODEL = "plearn_model"
KEY_PLEARN_HEADER = "plearn_header"
KEY_DESIGNATION_CACHE="designation_cache"
KEY_ASSIGNMENT_CACHE = "assignment_cache"

# Version KEYS
KEY_PANGOLEARN_VERSION="pangoLEARN_version"
KEY_PANGOLIN_VERSION="pangolin_version"
KEY_CONSTELLATIONS_VERSION="constellation_version"
KEY_SCORPIO_VERSION="scorpio_version"
KEY_PANGO_VERSION="pango_version"
KEY_PANGO_DESIGNATION_VERSION="pango_designation_version"

KEY_VERBOSE="verbose"
KEY_LOG_API = "log_api"
KEY_THREADS="threads"

# File names 
designation_cache_file = "lineages.hash.csv"

pangolearn_files = {
    "decisionTree_v1.joblib":KEY_PLEARN_MODEL,
    "decisionTreeHeaders_v1.joblib":KEY_PLEARN_HEADER
    }
usher_files = {
    "lineageTree.pb":KEY_USHER_PB
    }

# Dependencies
dependency_list = ["gofasta","minimap2","snakemake"]
module_list = ["Bio","sklearn","pandas","joblib","pysam","pangoLEARN","constellations"]
