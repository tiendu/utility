#!/bin/bash

# Help function
usage() {
cat << EOF
Usage:
    --input, -i
        absolute file path to fasta/fastq file e.g., $PWD/file.fasta or $PWD/file.fastq
        accept .gz files
    --output, -o (OPTIONAL)
        absolute file path to output
    --database, -d
        select either one of these: meta (Metagenomic), wgs (WGS), refseq (RefSeq), run (PH3HS_Runs), biosample (PH3HS_Biosample), hrna_2021 (Human RNAseq 2021)
Example:
    bash $0 -i $PWD/file.fasta -d refseq -o result.tsv
EOF
}

cmdline() {
    local arg=
    for arg; do
        local delim=""
        case "$arg" in
            # Translate --gnu-long-options to -g (short options)
            --input) args="${args}-i " ;;
            --output) args="${args}-o " ;;
            --database) args="${args}-d " ;;
            --help) args="${args}-h " ;;
            # Pass through anything else
            *) [[ "${arg:0:1}" == "-" ]] || delim="\""
                args="${args}${delim}${arg}${delim} " ;;
        esac
    done

    # Indexing option for later.
    OPTIND=1

    # Reset the positional parameters to the short options.
    eval set -- "$args"

    # Parse the options
    while getopts "i:o:d:h" OPTION; do
        case $OPTION in
            i) readonly INPUT_FILE="${OPTARG}" ;;
            o) readonly OUTPUT_FILE="${OPTARG}" ;;
            d) readonly DATABASE="${OPTARG}" ;;
            h) usage; exit 0 ;;
        esac
    done

    # Return help if no args.
    if [ $OPTIND -eq 1 ]; then
        usage
        exit 0
    fi

    # Check if both -i and -d are specified.
    if [ -z "${INPUT_FILE}" ] || [ -z "${DATABASE}" ]; then
        usage
        exit 1
    fi
}

pebblescout_request() {
    local sequence=$1
    local database=$2
    local nrow=$3
    local output=$4
    local g=10000
    local c=1000

    # Assemble the curl command into a variable
    local cmd="curl -s 'https://pebblescout.ncbi.nlm.nih.gov/sra-cl-be/sra-cl-be.cgi?rettype=pebblescout' \
        -H 'content-type: multipart/form-data; boundary=BOUNDARY' \
        --data-raw $'--BOUNDARY\\r\\nContent-Disposition: form-data; name=\"m\"\\r\\n\\r\\n2\\r\\n--BOUNDARY\\r\\nContent-Disposition: form-data; name=\"g\"\\r\\n\\r\\n${g}\\r\\n--BOUNDARY\\r\\nContent-Disposition: form-data; name=\"c\"\\r\\n\\r\\n${c}\\r\\n--BOUNDARY\\r\\nContent-Disposition: form-data; name=\"_r\"\\r\\n\\r\\n${nrow}\\r\\n--BOUNDARY\\r\\nContent-Disposition: form-data; name=\"accession\"\\r\\n\\r\\n\\r\\n--BOUNDARY\\r\\nContent-Disposition: form-data; name=\"_h\"\\r\\n\\r\\n1001\\r\\n--BOUNDARY\\r\\nContent-Disposition: form-data; name=\"fasta\"\\r\\n\\r\\n${sequence}\\r\\n--BOUNDARY\\r\\nContent-Disposition: form-data; name=\"from\"\\r\\n\\r\\n\\r\\n--BOUNDARY\\r\\nContent-Disposition: form-data; name=\"to\"\\r\\n\\r\\n\\r\\n--BOUNDARY\\r\\nContent-Disposition: form-data; name=\"db\"\\r\\n\\r\\n${database}\\r\\n--BOUNDARY\\r\\nContent-Disposition: form-data; name=\"retmode\"\\r\\n\\r\\n\\r\\n--BOUNDARY--\\r\\n' \
        --compressed"

    # Log the command
    echo "$cmd" | tee -a $LOGFILE

    # Execute the command
    eval "$cmd" 2>&1 | awk -v n="${nrow}" '/QueryID/ {for (i=0; i<n; i++) {getline; if ($0!~/^#.*$/ && $0!~/^$/) {print}}}' >> ${output}

    sleep 1
}

main() {
    cmdline "$@"

    # Check if file exists
    if [[ ! -f "${INPUT_FILE}" ]]; then
        echo "Error: ${INPUT_FILE} does not exist!"
        exit 1
    fi

    # Define output file
    if [ -z "${OUTPUT_FILE}" ]; then
        OUTPUT_FILE="${INPUT_FILE%.*}"
    fi

    echo "" > "${OUTPUT_FILE}"
    LOGFILE="${INPUT_FILE}_$(date +%Y%m%d%H%M%S).cmd"
    echo "" > $LOGFILE

    # Process input file
    if [[ "${INPUT_FILE}" =~ ".gz" ]]; then
        if [[ "${INPUT_FILE}" =~ ".fa" || "${INPUT_FILE}" =~ ".fasta" ]]; then
            zcat $INPUT_FILE | awk 'BEGIN {ORS=""} {if ($0 ~ /^>/) {printf "%s%s%s", (NR==1 ? "" : "\n"), $0, "\n"} else {print toupper($0)}}'  > ${INPUT_FILE}.tmp
        elif [[ "${INPUT_FILE}" =~ ".fq" || "${INPUT_FILE}" =~ ".fastq" ]]; then
            zcat $INPUT_FILE | awk 'NR%4==1 {sub(/^@/, ">", $0); print; getline; print}' > ${INPUT_FILE}.tmp
        else
            echo "Error: Unsupported file format!"
            exit 1
        fi
    else
        if [[ "${INPUT_FILE}" =~ ".fa" || "${INPUT_FILE}" =~ ".fasta" ]]; then
            awk 'BEGIN {ORS=""} {if ($0 ~ /^>/) {printf "%s%s%s", (NR==1 ? "" : "\n"), $0, "\n"} else {print toupper($0)}}' $INPUT_FILE > ${INPUT_FILE}.tmp
        elif [[ "${INPUT_FILE}" =~ ".fq" || "${INPUT_FILE}" =~ ".fastq" ]]; then
            awk 'NR%4==1 {sub(/^@/, ">", $0); print; getline; print}' $INPUT_FILE > ${INPUT_FILE}.tmp
        else
            echo "Error: Unsupported file format!"
            exit 1
        fi
    fi


    while read -r sequence; do
        echo ""
        pebblescout_request "$sequence" "${DATABASE}" 1 "${OUTPUT_FILE}"
    done < <(awk '/^>/ {gsub(/ /, "_", $0); printf $0 "\\r\\n"; getline; print}' ${INPUT_FILE}.tmp)

    # Clean up
    rm -f ${INPUT_FILE}.tmp
}

main "$@"
