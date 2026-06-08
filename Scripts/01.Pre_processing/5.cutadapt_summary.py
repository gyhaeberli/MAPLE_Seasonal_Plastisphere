#!/usr/bin/env python3
"""
Parse cutadapt log files and combine with Novogene raw counts
CORRECTED: Sum all replicates per sample before calculating retention
"""

import os
import re
import pandas as pd

# --- CONFIG ---
BASE_DIR = "/home/ahomew/g/gabrihae/data/glennsdata/MAPLE/18S/EUK_data/X204SC25074012-Z01-F001/01.RawData"
NOVOGENE_CSV = "/home/ahomew/g/gabrihae/data/glennsdata/MAPLE/18S/EUK_data/X204SC25074012-Z01-F001/NOVOGENE_raw_summary.csv"
LIBRARIES = ["GAB_EUK1", "GAB_EUK2", "GAB_EUK3", "GAB_EUK4"]
RUNS = ["run1", "run2"]
OUTPUT_FILE_DETAILED = "cutadapt_summary_detailed.csv"
OUTPUT_FILE_PER_SAMPLE = "cutadapt_summary_per_sample.csv"


def parse_cutadapt_log(log_file):
    """Parse a single cutadapt log file"""
    if not os.path.exists(log_file):
        return []
    
    results = []
    
    with open(log_file, 'r') as f:
        content = f.read()
    
    sample_sections = re.split(r'={5,}\s+(.+?)\s+={5,}', content)
    
    for i in range(1, len(sample_sections), 2):
        if i + 1 >= len(sample_sections):
            break
            
        sample_name = sample_sections[i].strip()
        section_content = sample_sections[i + 1]
        
        total_reads = None
        reads_written = None
        
        total_match = re.search(r'Total read pairs processed:\s+([\d,]+)', section_content)
        if total_match:
            total_reads = int(total_match.group(1).replace(',', ''))
        
        written_match = re.search(r'Pairs written \(passing filters\):\s+([\d,]+)', section_content)
        if written_match:
            reads_written = int(written_match.group(1).replace(',', ''))
        
        if total_reads is not None and reads_written is not None:
            results.append({
                'sample_basename': sample_name,
                'post_sabre_reads': total_reads,
                'post_cutadapt_reads': reads_written,
            })
    
    return results


def extract_sample_and_run_from_filename(sample_basename):
    """Extract sample name and run identifier"""
    mkdl_match = re.search(r'(MKDL\d+-\w+_\w+_L\d+)', sample_basename)
    if not mkdl_match:
        return None, None
    
    run_id = mkdl_match.group(1)
    sample_name = sample_basename.split(run_id)[0].rstrip('_')
    
    return sample_name, run_id


def main():
    print("="*60)
    print("Cutadapt Summary - CORRECTED for Replicate Summing")
    print("="*60)
    print()
    
    # Read Novogene raw counts
    print(f"Reading Novogene data from: {NOVOGENE_CSV}")
    novogene_df = pd.read_csv(NOVOGENE_CSV)
    novogene_df['Raw reads'] = novogene_df['Raw reads'] / 2  # Convert to pairs
    print("NOTE: Converted Novogene 'Raw reads' from individual reads to pairs")
    print(f"Novogene data: {len(novogene_df)} rows")
    
    all_data = []
    
    # Process each library and run
    for lib in LIBRARIES:
        for run in RUNS:
            work_dir = os.path.join(BASE_DIR, lib, run)
            
            if not os.path.isdir(work_dir):
                continue
            
            print(f"Processing {lib}/{run}...")
            
            for rep_num in range(1, 9):
                log_file = os.path.join(work_dir, f"rep{rep_num}_cutadapt.log")
                
                if not os.path.isfile(log_file):
                    continue
                
                sample_data = parse_cutadapt_log(log_file)
                
                for entry in sample_data:
                    sample_name, run_id = extract_sample_and_run_from_filename(entry['sample_basename'])
                    
                    if not sample_name or not run_id:
                        continue
                    
                    entry['library'] = lib
                    entry['run'] = run
                    entry['replicate'] = rep_num
                    entry['sample'] = sample_name
                    entry['run_identifier'] = run_id
                    all_data.append(entry)
    
    if not all_data:
        print("\nERROR: No data found!")
        return
    
    # Create detailed DataFrame (per replicate)
    detailed_df = pd.DataFrame(all_data)
    
    # Merge with Novogene for detailed view
    detailed_merged = detailed_df.merge(
        novogene_df,
        left_on=['sample', 'run_identifier'],
        right_on=['Sample', 'Library_Flowcell_Lane'],
        how='left'
    )
    
    # Save detailed (per-replicate) data
    detailed_merged.to_csv(os.path.join(BASE_DIR, OUTPUT_FILE_DETAILED), index=False)
    print(f"✓ Detailed data saved to: {OUTPUT_FILE_DETAILED}")
    
    # NOW: Sum all replicates for each sample
    print("\nSumming replicates per sample...")
    
    per_sample_df = detailed_df.groupby(['library', 'run', 'sample', 'run_identifier']).agg({
        'post_sabre_reads': 'sum',
        'post_cutadapt_reads': 'sum'
    }).reset_index()
    
    # Merge with Novogene
    merged_df = per_sample_df.merge(
        novogene_df,
        left_on=['sample', 'run_identifier'],
        right_on=['Sample', 'Library_Flowcell_Lane'],
        how='left'
    )
    
    # Calculate retention statistics
    merged_df['raw_reads'] = merged_df['Raw reads']
    merged_df['sabre_loss'] = merged_df['raw_reads'] - merged_df['post_sabre_reads']
    merged_df['cutadapt_loss'] = merged_df['post_sabre_reads'] - merged_df['post_cutadapt_reads']
    merged_df['total_loss'] = merged_df['raw_reads'] - merged_df['post_cutadapt_reads']
    
    merged_df['sabre_retention_%'] = (merged_df['post_sabre_reads'] / merged_df['raw_reads'] * 100).round(2)
    merged_df['cutadapt_retention_%'] = (merged_df['post_cutadapt_reads'] / merged_df['post_sabre_reads'] * 100).round(2)
    merged_df['total_retention_%'] = (merged_df['post_cutadapt_reads'] / merged_df['raw_reads'] * 100).round(2)
    
    # Save per-sample summary
    output_cols = [
        'library', 'run', 'sample',
        'raw_reads', 'post_sabre_reads', 'post_cutadapt_reads',
        'sabre_loss', 'cutadapt_loss', 'total_loss',
        'sabre_retention_%', 'cutadapt_retention_%', 'total_retention_%'
    ]
    
    final_df = merged_df[output_cols].copy()
    final_df = final_df.sort_values(['library', 'run', 'sample'])
    
    output_path = os.path.join(BASE_DIR, OUTPUT_FILE_PER_SAMPLE)
    final_df.to_csv(output_path, index=False)
    print(f"✓ Per-sample summary saved to: {OUTPUT_FILE_PER_SAMPLE}")
    
    # Print summary statistics
    print("\n" + "="*60)
    print("SUMMARY STATISTICS (PER SAMPLE, ALL REPS SUMMED)")
    print("="*60)
    
    total_raw = final_df['raw_reads'].sum()
    total_post_sabre = final_df['post_sabre_reads'].sum()
    total_post_cutadapt = final_df['post_cutadapt_reads'].sum()
    
    print(f"\nTotal samples: {len(final_df)}")
    print(f"Total raw reads (Novogene): {total_raw:,}")
    print(f"Total post-sabre reads (all 8 reps): {total_post_sabre:,}")
    print(f"Total post-cutadapt reads: {total_post_cutadapt:,}")
    print(f"\nSabre retention: {total_post_sabre/total_raw*100:.2f}%")
    print(f"Cutadapt retention: {total_post_cutadapt/total_post_sabre*100:.2f}%")
    print(f"TOTAL retention (raw → final): {total_post_cutadapt/total_raw*100:.2f}%")
    
    # By library
    print("\n" + "-"*60)
    print("BY LIBRARY:")
    print("-"*60)
    lib_summary = final_df.groupby('library').agg({
        'raw_reads': 'sum',
        'post_sabre_reads': 'sum',
        'post_cutadapt_reads': 'sum'
    })
    lib_summary['total_retention_%'] = (lib_summary['post_cutadapt_reads'] / lib_summary['raw_reads'] * 100).round(2)
    print(lib_summary.to_string())
    
    print("\n" + "="*60)
    print("Analysis complete!")
    print("="*60)


if __name__ == "__main__":
    main()