import os
import sys
import subprocess

def install_dependencies():
    """Install required Python packages"""
    try:
        import Bio
        print("✅ BioPython already installed")
    except ImportError:
        print("📦 Installing BioPython...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython"])
        print("✅ BioPython installed")

# Install dependencies before importing Bio modules
install_dependencies()

from Bio import SeqIO

# ==================== KAGGLE CONFIGURATION ====================
# This script compares R DADA2 vs Python DADA2 outputs
# Validates that both implementations produce similar results

dada2_output_path = "/kaggle/input/datasets/aramelheni/16s-rrna-microbiome-sequencing-analysis-dataset/filtered/results/dada2_output.fastq"
python_output_path = "/kaggle/input/datasets/aramelheni/16s-rrna-microbiome-sequencing-analysis-dataset/filtered/results/python_output.fastq"

print("=" * 60)
print("DADA2 Output Comparison: R vs Python")
print("=" * 60)
print(f"📁 R DADA2 output: {dada2_output_path}")
print(f"📁 Python output: {python_output_path}")
print()

# Check if files exist
if not os.path.exists(dada2_output_path):
    print(f"❌ Error: R DADA2 output not found: {dada2_output_path}")
    exit(1)
if not os.path.exists(python_output_path):
    print(f"❌ Error: Python output not found: {python_output_path}")
    exit(1)

# Function to read FASTQ into list of (id, sequence, quality)
def read_fastq(file_path):
    records = []
    skipped = 0
    line_num = 0
    
    with open(file_path, "r") as handle:
        while True:
            try:
                # FASTQ format: 4 lines per record
                # Line 1: @ID
                # Line 2: Sequence
                # Line 3: +
                # Line 4: Quality scores
                
                header = handle.readline()
                line_num += 1
                if not header:
                    break  # End of file
                
                if not header.startswith('@'):
                    print(f"⚠️  Skipping malformed record at line {line_num} in {os.path.basename(file_path)}")
                    skipped += 1
                    continue
                
                seq = handle.readline().strip()
                line_num += 1
                plus = handle.readline()
                line_num += 1
                qual = handle.readline().strip()
                line_num += 1
                
                # Validate lengths match
                if len(seq) != len(qual):
                    print(f"⚠️  Skipping record {header.strip()} - sequence length ({len(seq)}) != quality length ({len(qual)})")
                    skipped += 1
                    continue
                
                # Convert quality string to Phred scores
                qual_scores = [ord(c) - 33 for c in qual]  # Assuming Phred+33 encoding
                
                record_id = header.strip()[1:]  # Remove @ prefix
                records.append((record_id, seq, qual_scores))
                
            except Exception as e:
                print(f"⚠️  Error reading record at line {line_num} in {os.path.basename(file_path)}: {str(e)}")
                skipped += 1
                continue
    
    if skipped > 0:
        print(f"   Total skipped: {skipped} problematic record(s)")
    
    return records

# Load both files
print("📖 Reading output files...\n")
dada2_records = read_fastq(dada2_output_path)
python_records = read_fastq(python_output_path)

# Compare number of sequences
print("📊 Sequence Statistics:")
print(f"   R DADA2 ASVs: {len(dada2_records)}")
print(f"   Python ASVs: {len(python_records)}")

if len(dada2_records) != len(python_records):
    print(f"\n⚠️  Different number of ASVs detected")
    print(f"   Difference: {abs(len(dada2_records) - len(python_records))} sequences\n")
else:
    print(f"\n✅ Same number of ASVs: {len(dada2_records)}\n")

# Compare sequences
print("🔍 Comparing sequences...\n")

# Create sets of sequences for comparison
dada2_seqs = set(seq for _, seq, _ in dada2_records)
python_seqs = set(seq for _, seq, _ in python_records)

# Find common and unique sequences
common_seqs = dada2_seqs & python_seqs
only_in_dada2 = dada2_seqs - python_seqs
only_in_python = python_seqs - dada2_seqs

print(f"✅ Common sequences: {len(common_seqs)}")
print(f"📍 Only in R DADA2: {len(only_in_dada2)}")
print(f"📍 Only in Python: {len(only_in_python)}")

if len(common_seqs) == len(dada2_seqs) == len(python_seqs):
    print("\n🎉 Perfect match! All sequences are identical between R and Python implementations!")
elif len(common_seqs) > 0:
    similarity = (len(common_seqs) / max(len(dada2_seqs), len(python_seqs))) * 100
    print(f"\n📈 Sequence similarity: {similarity:.1f}%")
    print("\n⚠️  Note: Some differences are expected due to:")
    print("   - Different denoising algorithms (R DADA2 is more sophisticated)")
    print("   - Different error model implementations")
    print("   - Random seed differences in clustering")
else:
    print("\n❌ No common sequences found. This suggests a major implementation difference.")

# Compare quality scores for common sequences
if len(common_seqs) > 0:
    print("\n📊 Quality Score Comparison:")
    dada2_dict = {seq: qual for _, seq, qual in dada2_records}
    python_dict = {seq: qual for _, seq, qual in python_records}
    
    qual_matches = 0
    total_compared = 0
    for seq in list(common_seqs)[:5]:  # Check first 5 common sequences
        if seq in dada2_dict and seq in python_dict:
            dada2_qual = dada2_dict[seq]
            python_qual = python_dict[seq]
            if dada2_qual == python_qual:
                qual_matches += 1
            total_compared += 1
    
    if total_compared > 0:
        print(f"   Quality scores match: {qual_matches}/{total_compared} sequences checked")
        if qual_matches < total_compared:
            print("   ⚠️  Quality scores may differ due to different merging algorithms")

print("\n" + "=" * 60)
print("Comparison complete!")
print("=" * 60)
