## Novogene ran 2 separate sequencing runs on each samples. We will treat them as different entities and merge them later. 
## This script will organize the files according to their sequencing run.

#!/usr/bin/env python3
import os
import shutil

# --- CONFIG ---
BASE_DIR = "/home/ahomew/g/gabrihae/data/glennsdata/MAPLE/18S/EUK_data/X204SC25074012-Z01-F001/01.RawData"

# Library-specific run identifiers
LIBRARY_CONFIG = {
    "GAB_EUK1": {
        "run1": "MKDL250008974-1A_22WW55LT4",
        "run2": "MKDL250008974-1A_22YHHNLT4"
    },
    "GAB_EUK2": {
        "run1": "MKDL250008975-1A_22WW55LT4",
        "run2": "MKDL250008975-1A_22YJLGLT4"
    },
    "GAB_EUK3": {
        "run1": "MKDL250008976-1A_22WW55LT4",
        "run2": "MKDL250008976-1A_22YJLGLT4"
    },
    "GAB_EUK4": {
        "run1": "MKDL250008977-1A_22WW55LT4",
        "run2": "MKDL250008977-1A_22YJLGLT4"
    }
}

ACCEPTED_EXTS = (".fq.gz", ".fastq.gz")
DO_MOVE = True   # True: move files; False: create symlinks


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def is_fastq(name):
    return name.endswith(ACCEPTED_EXTS)


def process_library(gab_dir, lib_name, run1_token, run2_token):
    """Process a single library folder"""
    print(f"\n{'='*60}")
    print(f"Processing library: {lib_name}")
    print(f"  Run1 token: {run1_token}")
    print(f"  Run2 token: {run2_token}")
    print(f"{'='*60}")
    
    if not os.path.isdir(gab_dir):
        print(f"WARNING: Directory not found: {gab_dir}")
        print("Skipping this library...")
        return None
    
    run1_dir = os.path.join(gab_dir, "run1")
    run2_dir = os.path.join(gab_dir, "run2")
    ensure_dir(run1_dir)
    ensure_dir(run2_dir)
    
    moved_run1 = moved_run2 = skipped = 0
    
    # Look only one level down (sample folders)
    for entry in sorted(os.listdir(gab_dir)):
        sample_dir = os.path.join(gab_dir, entry)
        if not os.path.isdir(sample_dir):
            continue
        if entry in {"run1", "run2"}:
            continue  # don't recurse into the run folders themselves
        
        # Scan files in this sample folder
        for fname in sorted(os.listdir(sample_dir)):
            if not is_fastq(fname):
                continue
            
            src = os.path.join(sample_dir, fname)
            
            if run1_token in fname:
                dest = os.path.join(run1_dir, fname)
                if DO_MOVE:
                    shutil.move(src, dest)
                else:
                    if not os.path.exists(dest):
                        os.symlink(os.path.relpath(src, run1_dir), dest)
                moved_run1 += 1
            
            elif run2_token in fname:
                dest = os.path.join(run2_dir, fname)
                if DO_MOVE:
                    shutil.move(src, dest)
                else:
                    if not os.path.exists(dest):
                        os.symlink(os.path.relpath(src, run2_dir), dest)
                moved_run2 += 1
            
            else:
                skipped += 1
    
    print(f"\nResults for {lib_name}:")
    print(f"  Moved to run1: {moved_run1}")
    print(f"  Moved to run2: {moved_run2}")
    print(f"  Skipped (no token match): {skipped}")
    
    return moved_run1, moved_run2, skipped


def main():
    if not os.path.isdir(BASE_DIR):
        raise SystemExit(f"Base directory not found: {BASE_DIR}")
    
    total_run1 = total_run2 = total_skipped = 0
    
    # Process each library with its specific run tokens
    for lib_name, tokens in LIBRARY_CONFIG.items():
        lib_path = os.path.join(BASE_DIR, lib_name)
        results = process_library(
            lib_path, 
            lib_name, 
            tokens["run1"], 
            tokens["run2"]
        )
        if results:
            total_run1 += results[0]
            total_run2 += results[1]
            total_skipped += results[2]
    
    # Summary
    print(f"\n{'='*60}")
    print("OVERALL SUMMARY")
    print(f"{'='*60}")
    print(f"Total files moved to run1: {total_run1}")
    print(f"Total files moved to run2: {total_run2}")
    print(f"Total files skipped: {total_skipped}")
    print(f"\nAll {len(LIBRARY_CONFIG)} libraries processed!")


if __name__ == "__main__":
    main()