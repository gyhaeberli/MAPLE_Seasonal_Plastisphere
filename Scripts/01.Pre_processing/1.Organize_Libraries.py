# This script is an organizer: moves sequencing sample folders into specific library subfolders (EUK1 samples will go into a GAB_EUK1 folder)
# Here we need a CSV instruction file -- assigning samples to libraries (see Novogene order sheet)


import os
import pandas as pd
import shutil

# --- CONFIG: File paths ---
CSV_PATH = "/home/ahomew/g/gabrihae/data/glennsdata/MAPLE/18S/EUK_data/X204SC25074012-Z01-F001/Processing_scripts/EUK_Library_sorting.csv"
SOURCE_ROOT = "/home/ahomew/g/gabrihae/data/glennsdata/MAPLE/18S/EUK_data/X204SC25074012-Z01-F001/01.RawData"


# --- LOAD & CLEAN ---
# Read all columns as string and don't auto-convert blanks/NA-like text to NaN
df = pd.read_csv(CSV_PATH, dtype=str, keep_default_na=False)

# Ensure required columns exist (adjust names here if your CSV uses different headers)
required_cols = {"sample", "library"}
missing_cols = required_cols - set(df.columns.str.lower())
if missing_cols:
    raise ValueError(f"CSV missing required columns: {missing_cols} (have: {list(df.columns)})")

# Normalize to exact column names
# If your CSV headers are 'Sample'/'Library', this maps them
cols_lower = {c.lower(): c for c in df.columns}
sample_col = cols_lower["sample"]
library_col = cols_lower["library"]

# Strip whitespace
df[sample_col]  = df[sample_col].astype(str).str.strip()
df[library_col] = df[library_col].astype(str).str.strip()

# Drop empty rows
df = df[(df[sample_col] != "") & (df[library_col] != "")].copy()

# Remove any path separators just in case
df[sample_col]  = df[sample_col].str.replace(r"[\\/]+", "_", regex=True)
df[library_col] = df[library_col].str.replace(r"[\\/]+", "_", regex=True)

# --- ENSURE DESTINATION FOLDERS ---
libraries = sorted(df[library_col].unique())
for lib in libraries:
    dest_path = os.path.join(SOURCE_ROOT, lib)
    os.makedirs(dest_path, exist_ok=True)

# --- MOVE ---
missing, moved, skipped = [], [], []

for _, row in df.iterrows():
    sample = row[sample_col]
    library = row[library_col]

    sample_dir = os.path.join(SOURCE_ROOT, sample)
    dest_dir   = os.path.join(SOURCE_ROOT, library, sample)

    if not os.path.exists(sample_dir):
        missing.append(sample)
        print(f"Missing: {sample}")
        continue

    # If already moved previously, don't crash
    if os.path.exists(dest_dir):
        skipped.append(sample)
        print(f"Already in place: {sample}")
        continue

    shutil.move(sample_dir, dest_dir)
    moved.append(sample)
    print(f"Moved {sample} → {library}")

# --- SUMMARY ---
print("\n--- Summary ---")
print(f"Moved: {len(moved)}")
print(f"Already in place: {len(skipped)}")
print(f"Missing: {len(missing)}")
if missing:
    print("Examples of missing:", ", ".join(missing[:10]))
