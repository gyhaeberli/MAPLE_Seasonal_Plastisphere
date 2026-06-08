## Because we had 2 sequencing runs for each library, we divided these runs in separate directories: run1 and run2 
# this was done using EUK_orgRuns.py

##This bash script will demultiplex the sequencing data that is in one directory. Seq files need to be in the same place.

##Ultimately, this script will write another .sh script to demultiplex the data for each run.
###########################################################################################################
#Change to the directory containing the sequencing data of 1 run at the time
#cd GAB_EUK1/run1
#cd GAB_EUK1/run2

#cd GAB_EUK2/run1
#cd GAB_EUK2/run2

#cd GAB_EUK3/run1
#cd GAB_EUK3/run2

#cd GAB_EUK4/run1
cd GAB_EUK4/run2


## Seq files all in one place - now demtuliplexing using sabre
echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz}; 
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /data/glennsdata/MAPLE/18S/EUK_data/X204SC25074012-Z01-F001/Processing_scripts/barcode_18_EUKforward.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
mv rep1f ${bn}_rep1f.fq;
mv rep1r ${bn}_rep1r.fq;
mv rep2f ${bn}_rep2f.fq;
mv rep2r ${bn}_rep2r.fq;
mv rep3f ${bn}_rep3f.fq;
mv rep3r ${bn}_rep3r.fq;
mv rep4f ${bn}_rep4f.fq;
mv rep4r ${bn}_rep4r.fq
mv rep1f_2 ${bn}_rep1f_2.fq;
mv rep1r_2 ${bn}_rep1r_2.fq;
mv rep2f_2 ${bn}_rep2f_2.fq;
mv rep2r_2 ${bn}_rep2r_2.fq;
mv rep3f_2 ${bn}_rep3f_2.fq;
mv rep3r_2 ${bn}_rep3r_2.fq;
mv rep4f_2 ${bn}_rep4f_2.fq;
mv rep4r_2 ${bn}_rep4r_2.fq;done' > maple18s_EUK1.sh

bash maple18s_EUK1.sh

## copy paste this in terminal!


## Check if it worked properly

cd /home/ahomew/g/gabrihae/data/glennsdata/MAPLE/18S/EUK_data/X204SC25074012-Z01-F001/01.RawData

for lib in GAB_EUK1 GAB_EUK2 GAB_EUK3 GAB_EUK4; do
    for run in run1 run2; do
        if [ -d "$lib/$run" ]; then
            echo "=== $lib/$run ==="
            cd "$lib/$run"
            
            input=$(($(zcat *_1.fq.gz 2>/dev/null | wc -l) / 4))
            output=$(($(cat *_rep*f.fq *_rep*f_2.fq *unassigned*.fq 2>/dev/null | wc -l) / 4))
            
            if [ $input -gt 0 ]; then
                retention=$(awk "BEGIN {printf \"%.2f\", ($output/$input)*100}")
                echo "  Input:  $input pairs"
                echo "  Output: $output pairs"
                echo "  Retention: $retention%"
                
                if (( $(echo "$retention < 95" | bc -l) )); then
                    echo "  ⚠ WARNING: Low retention!"
                elif (( $(echo "$retention > 105" | bc -l) )); then
                    echo "  ⚠ WARNING: High retention (check for double-counting)"
                else
                    echo "  ✓ Good"
                fi
            else
                echo "  (No data)"
            fi
            
            cd ../..
            echo ""
        fi
    done
done

####### LAST OUTPUT ######

     === GAB_EUK1/run1 ===
    Input:  11215276 pairs
    Output: 11379500 pairs
    Retention: 101.46%
    ✓ Good

    === GAB_EUK1/run2 ===
    Input:  26435805 pairs
    Output: 26807575 pairs
    Retention: 101.41%
    ✓ Good

    === GAB_EUK2/run1 ===
    Input:  11474651 pairs
    Output: 11678279 pairs
    Retention: 101.77%
    ✓ Good

    === GAB_EUK2/run2 ===
    Input:  23425472 pairs
    Output: 23754954 pairs
    Retention: 101.41%
    ✓ Good

    === GAB_EUK3/run1 ===
    Input:  4372030 pairs
    Output: 4430507 pairs
    Retention: 101.34%
    ✓ Good

    === GAB_EUK3/run2 ===
    Input:  24011531 pairs
    Output: 24353525 pairs
    Retention: 101.42%
    ✓ Good

    === GAB_EUK4/run1 ===
    Input:  5581633 pairs
    Output: 5658091 pairs
    Retention: 101.37%
    ✓ Good

    === GAB_EUK4/run2 ===
    Input:  15168728 pairs
    Output: 15390842 pairs
    Retention: 101.46%
    ✓ Good
    
