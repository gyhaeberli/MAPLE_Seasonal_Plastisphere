### 
#### CUTADAPT

## now using cutadapt to find the RV internal index for that PCR rep 
# and fully dereplicate, while also removing any chimeras 
##base eukaryotic primers are:

# 1391Feu_F (fwd): GTACACACCGCCCGTC
# EukBeu_R (rcRv): GTAGGTGAACCTGCAGAAGGATCA
# 1391Feu_F (rcfwd): GACGGGCGGTGTGTAC
# EukBeu_R (Rv): TGATCCTTCTGCAGGTTCACCTAC

#1	AGGAA
#2	GAGTGG
#3	CCACGTC
#4	TTCTCAGC
#5	CTAGG
#6	GCTTAT
#7	GCGAAGT
#8	AATCCTAT



#-a fwd...rcRv
# -A Rv...rcfwd

#!/bin/bash
set -euo pipefail

###########################################
# SET THIS TO WHERE YOUR FASTQ FILES ARE
###########################################
DIR="/home/ahomew/g/gabrihae/data/glennsdata/MAPLE/18S/EUK_data/X204SC25074012-Z01-F001/01.RawData/GAB_EUK4/run2"
# change directory for all libraries/run --> multiple screens

cd "$DIR" || {
    echo "ERROR: Could not cd into directory:"
    echo "  $DIR"
    exit 1
}

echo "Working directory: $DIR"
echo

###########################################
# CONSTANT PRIMERS
###########################################

# Forward primer (5' of R1)
FWD="GTACACACCGCCCGTC"

# Reverse primer (as seen in R2, without barcode)
RV_BASE="TGATCCTTCTGCAGGTTCACCTAC"

# Reverse complement of forward primer (internal trimming on R1/R2)
RC_FWD="GACGGGCGGTGTGTAC"

###########################################
# FUNCTION TO PROCESS ONE REPLICATE
###########################################
# run_rep <rep_number> <R1_pattern> <R2_pattern> <R2_primer_with_barcode>
run_rep () {
    local rep="$1"
    local R1pat="$2"
    local R2pat="$3"
    local R2primer_bc="$4"

    local log="rep${rep}_cutadapt.log"
    : > "$log"

    echo "=== Processing Replicate ${rep} ==="

    for f in *"$R1pat"; do
        [[ -e "$f" ]] || continue

        # Derive basename by stripping the R1 pattern
        bn="${f/$R1pat/}"
        r="${bn}${R2pat}"

        if [[ ! -f "$r" ]]; then
            echo "WARNING: Missing R2 for $f (expected $r)" | tee -a "$log"
            continue
        fi

        echo "Cutadapt (rep${rep}): R1=$f  R2=$r" | tee -a "$log"

        # Temporary file to capture this run's cutadapt output
        tmp=$(mktemp)

        # Logic:
        #  -g FWD           : R1 must contain forward primer at 5'
        #  -G R2primer_bc   : R2 must contain barcoded reverse primer at 5'
        #  -a RC_FWD        : trim any occurrence of rev-comp forward primer internally/3'
        #  -A RV_BASE       : trim any occurrence of base reverse primer internally/3'
        #  --pair-filter=both + --discard-untrimmed:
        #       keep only pairs where BOTH reads had their assigned primer
        cutadapt \
            -g "$FWD" \
            -G "$R2primer_bc" \
            -a "$RC_FWD" \
            -A "$RV_BASE" \
            --pair-filter=both \
            --discard-untrimmed \
            -o "${bn}.rep${rep}.trim1.fq.gz" \
            -p "${bn}.rep${rep}.trim2.fq.gz" \
            "$f" "$r" \
            > "$tmp" 2>&1

        {
            echo "===== ${bn} ====="
            cat "$tmp"
            echo
        } >> "$log"

        rm "$tmp"
    done

    echo "Replicate ${rep} complete."
    echo
}

###########################################
# RUN ALL 8 REPLICATES
###########################################

# From your notes:
# 1391Feu_F (fwd):   GTACACACCGCCCGTC
# 1391Feu_F (rcfwd): GACGGGCGGTGTGTAC
# EukBeu_R (Rv):     TGATCCTTCTGCAGGTTCACCTAC
# EukBeu_R (rcRv):   GTAGGTGAACCTGCAGAAGGATCA

######## Replicate 1 ########
# rep1 reverse barcode: AGGAA
# Rv+barcode (R2): AGGAATGATCCTTCTGCAGGTTCACCTAC
run_rep 1 "rep1f.fq" "rep1r.fq" "AGGAATGATCCTTCTGCAGGTTCACCTAC"

######## Replicate 2 ########
# rep2 reverse barcode: GAGTGG
# Rv+barcode (R2): GAGTGGTGATCCTTCTGCAGGTTCACCTAC
run_rep 2 "rep2f.fq" "rep2r.fq" "GAGTGGTGATCCTTCTGCAGGTTCACCTAC"

######## Replicate 3 ########
# rep3 reverse barcode: CCACGTC
# Rv+barcode (R2): CCACGTCTGATCCTTCTGCAGGTTCACCTAC
run_rep 3 "rep3f.fq" "rep3r.fq" "CCACGTCTGATCCTTCTGCAGGTTCACCTAC"

######## Replicate 4 ########
# rep4 reverse barcode: TTCTCAGC
# Rv+barcode (R2): TTCTCAGCTGATCCTTCTGCAGGTTCACCTAC
run_rep 4 "rep4f.fq" "rep4r.fq" "TTCTCAGCTGATCCTTCTGCAGGTTCACCTAC"

######## Replicate 5 ########
# rep5 reverse barcode: CTAGG
# Rv+barcode (R2): CTAGGTGATCCTTCTGCAGGTTCACCTAC
# Filenames: *rep1f_2.fq / *rep1r_2.fq
run_rep 5 "rep1f_2.fq" "rep1r_2.fq" "CTAGGTGATCCTTCTGCAGGTTCACCTAC"

######## Replicate 6 ########
# rep6 reverse barcode: GCTTAT
# Rv+barcode (R2): GCTTATTGATCCTTCTGCAGGTTCACCTAC
run_rep 6 "rep2f_2.fq" "rep2r_2.fq" "GCTTATTGATCCTTCTGCAGGTTCACCTAC"

######## Replicate 7 ########
# rep7 reverse barcode: GCGAAGT
# Rv+barcode (R2): GCGAAGTTGATCCTTCTGCAGGTTCACCTAC
run_rep 7 "rep3f_2.fq" "rep3r_2.fq" "GCGAAGTTGATCCTTCTGCAGGTTCACCTAC"

######## Replicate 8 ########
# rep8 reverse barcode: AATCCTAT
# Rv+barcode (R2): AATCCTATTGATCCTTCTGCAGGTTCACCTAC
run_rep 8 "rep4f_2.fq" "rep4r_2.fq" "AATCCTATTGATCCTTCTGCAGGTTCACCTAC"

echo "All replicates finished."






































