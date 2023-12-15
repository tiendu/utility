#!/bin/bash

# Help section here.
usage() {
cat << EOF
    Usage:
        --input, -i
            absolute file path e.g., $PWD/
        --extension, -e
            file extension e.g., "_R[1,2].fastq.gz"
        --primers, -p
            forward primer and reverse primer sequence e.g., "ATN GNN"
    Example:
        bash $0 -i $PWD/ -e "*_R[1,2].fastq.gz" -p "ATGN ACYN"
EOF
}

# Define input parameters.
cmdline() {
    local arg=
    for arg; do
        local delim=""
        case "$arg" in
            # Translate --gnu-long-options to -g (short options)
            --input) args="${args}-i " ;;
            --extension) args="${args}-e " ;;
            --primers) args="${args}-p " ;;
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

    while getopts "i:e:p:h" OPTION; do
        case $OPTION in
            i) readonly INPUT_PATH="${OPTARG}" ;;
            e) readonly FILE_EXTENSION="${OPTARG}" ;;
            p) readonly PRIMERS="${OPTARG}" ;;
            h) usage; exit 0 ;;
        esac
    done

    # Return help if no args and if args not supplemented enough.
    shift "$((OPTIND-1))"

    if [ $OPTIND -eq 1 ] || [ $OPTIND -lt 7 ]; then
        usage
        exit 0
    fi
}

check_format_fastq() {
    local line=''
    local file=$1
    while mapfile -t -n 4 line && ((${#line[@]})); do
        # By order, unless 1st line starts with @, 2nd line is comprised of ATGC, 3rd line starts with '+', 4th line is not empty, report the fastq file, and length of 2nd line must equal length of 4th line.
        if [[ ${line[0]} == @* ]] && [[ ${line[1]} == +(A|T|G|C|N) ]] && [[ ${line[2]} == \+* ]] && ! [ -z "${line[3]}" ] && [ ${#line[1]} -eq ${#line[3]} ]; then
            continue
        else
            echo "$file is not properly constructed!" >&2
            exit 1
        fi
    done < <(cat "$file")
}

delineate_sequence() {
    local sequence=$1
    awk '
    function delineate(s) {
        hash["A"]="A";
        hash["T"]="T";
        hash["G"]="G";
        hash["C"]="C";
        hash["R"]="A G";
        hash["Y"]="C T";
        hash["S"]="G C";
        hash["W"]="A T";
        hash["K"]="G T";
        hash["M"]="A C";
        hash["B"]="C G T";
        hash["D"]="A G T";
        hash["H"]="A C T";
        hash["V"]="A C G";
        hash["N"]="A T G C";
        solutions=1;
        tmp=s;
        m=split(tmp, a, "");
        for (i=1; i<=m; i++) {
            n=split(hash[a[i]], b, " ");
            solutions*=n;
            for (j=0; j<n; j++) {
                list[i][j]=b[j+1]
                }
            };
        for (i=1; i<=solutions; i++) {
            tmp="";
            idx=i;
            for (j=1; j<=m; j++) {
                tmp=tmp list[j][idx%length(list[j])];
                idx=int(idx/length(list[j]))
                };
            print tmp
            }
        } {delineate(toupper($0))}' <(echo $sequence)
}

reverse_complement() {
    local sequence=$1
    printf "%s" $sequence | tr "ATGC" "TACG" | rev
}

primer_check() {
    local sequence=$1
    local primer=$2
    awk -v FS=' ' '{
        sequence=$1;
        primer=$2;
        match(sequence, primer);
        if (RSTART==1) {
            print "FORWARD"
        } else if (RSTART+RLENGTH-1==length(sequence)){
            print "REVERSE"
        } else {
            print "NONE"
        }
    }' <(echo $sequence $primer)
}

main() {
    cmdline "$@"
    declare -A readonly names

    # Check if path exists and is properly formatted.
    if [ ! -d "${INPUT_PATH}" ]; then
        echo "${INPUT_PATH} not exists!" >&2
        exit 1
    fi
    if [[ ! "${INPUT_PATH}" == *\/ ]]; then
        echo "Missing / at the end of input path!" >&2
        exit 1
    fi

    # Check if fastq(s) exists.
    if ! compgen -G "${INPUT_PATH}*${FILE_EXTENSION}" > /dev/null; then
        echo "${INPUT_PATH}*${FILE_EXTENSION} not exists!" >&2
        exit 1
    fi


    # Check whether fastq(s) is properly formatted.
    for file in $(find "${INPUT_PATH}" -maxdepth 1 -name "*${FILE_EXTENSION}"); do
        check_format_fastq "$file"
        local name=$(echo "$file" | rev | cut -d'/' -f1 | rev | sed "s/${FILE_EXTENSION}//g")
        ((names[$name]=names[$name]+1))
    done

    # Check if paired-ends.
    for index in "${!names[@]}"; do
        if [ "${names[$index]}" -eq 2 ]; then
            names[$index]="PE"
        elif [ "${names[$index]}" -eq 1 ]; then
            names[$index]="SE"
        fi
    done

    # Delineate primer sequences.
    local forward_primer=$(echo $PRIMERS | cut -d' ' -f1)
    local reverse_primer=$(echo $PRIMERS | cut -d' ' -f2)
    forward_primer=$(delineate_sequence $forward_primer)
    reverse_primer=$(delineate_sequence $reverse_primer)

    # Reverse complement of reverse primer.
    local copy
    for i in ${reverse_primer[@]}; do
        local temporary=$(reverse_complement $i)
        copy=(${copy[@]} $temporary)
    done
    reverse_primer=${copy[@]}
    unset copy

    # Checking if the primers match.
    for index in "${!names[@]}"; do
        # Counting number of sequences and number of matches.
        local count_sequence=0
        local count_match=0
        for file in $(find "${INPUT_PATH}" -maxdepth 1 -name "${index}*${FILE_EXTENSION}"); do
            for sequence in `awk '/^@/ {getline; print}' $file`; do
                (( count_sequence++ ))
                for primer in ${reverse_primer[@]}; do
                    echo $primer
                    local check=$(primer_check $sequence $primer)
                    if [ "${check}" == "REVERSE" ]; then
                        (( count_match++ ))
                        break
                    fi
                done
                for primer in ${forward_primer[@]}; do
                    local check=$(primer_check $sequence $primer)
                    if [ "${check}" == "FORWARD" ]; then
                        (( count_match++ ))
                        break
                    fi
                done
            done
        done

        # For paired-end, the number of matches is doubled the number of sequences.
        if [ "${names[$index]}" == "PE" ]; then
            if [ ! $count_sequence -eq $(( $count_match / 2 )) ]; then
                echo "Primers don't match with ${index}!"
                unset 'names[$index]'
                exit 1
            fi
        # For single-end, the number of matches equals the number of sequences.
        elif [ "${names[$index]}" == "SE" ]; then
            if [ $count_sequence -eq $count_match ]; then
                echo "Primers don't match with ${index}!"
                unset 'names[$index]'
                exit 1
            fi
        fi
    done
}

main "$@"
