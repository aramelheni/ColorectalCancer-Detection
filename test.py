# %% [code]
import os
import subprocess
import sys

# ==================== KAGGLE CONFIGURATION ====================
# Download SRR29923448 FASTQ files from NCBI SRA and upload to Kaggle dataset
forward_file = "/kaggle/input/datasets/aramelheni/crc-microbiome-raw-fastq/SRR29923448_1.fastq"
reverse_file = "/kaggle/input/datasets/aramelheni/crc-microbiome-raw-fastq/SRR29923448_2.fastq"

# Download SILVA v138.1 from https://zenodo.org/record/4587955
# Extract .gz files and upload to Kaggle dataset
silva_folder = "/kaggle/input/datasets/aramelheni/silva-v138-reference"

# Output folders in Kaggle working directory
output_base = "/kaggle/working"
project_folder = output_base
r_script = "/kaggle/usr/lib/notebooks/aramelheni/dada2-preprocessing-pipeline-r/dada2_preprocessing_pipeline_r.R"
rscript_path = "Rscript"

# Create output folders
filtered_folder = os.path.join(output_base, "filtered")
reads_folder = os.path.join(filtered_folder, "reads")
plots_folder = os.path.join(filtered_folder, "plots")
results_folder = os.path.join(filtered_folder, "results")

for folder in [filtered_folder, reads_folder, plots_folder, results_folder]:
    os.makedirs(folder, exist_ok=True)
    
print(f"📁 Forward reads: {forward_file}")
print(f"📁 Reverse reads: {reverse_file}")
print(f"📁 SILVA database: {silva_folder}")
print(f"📁 Output directory: {output_base}")
print()

# Verify required files exist
if not os.path.exists(forward_file):
    print(f"❌ Error: Forward file not found: {forward_file}")
    sys.exit(1)
if not os.path.exists(reverse_file):
    print(f"❌ Error: Reverse file not found: {reverse_file}")
    sys.exit(1)
if not os.path.exists(silva_folder):
    print(f"❌ Error: SILVA folder not found: {silva_folder}")
    sys.exit(1)
    
print("🚀 Starting DADA2 pipelines...\n")

# Run R DADA2 script
print("=" * 60)
print("Running R DADA2 Pipeline")
print("=" * 60)
try:
    subprocess.run([
        rscript_path, r_script,
        forward_file, reverse_file,
        reads_folder, plots_folder, results_folder, silva_folder
    ], check=True)
    print("\n✅ R DADA2 pipeline complete!\n")
except subprocess.CalledProcessError as e:
    print(f"\n❌ R pipeline failed with error code {e.returncode}")
    sys.exit(1)

# Run Python DADA2 script
print("=" * 60)
print("Running Python DADA2 Pipeline")
print("=" * 60)
python_script = "/kaggle/usr/lib/notebooks/aramelheni/dada2-preprocessing-pipeline-python/dada2_preprocessing_pipeline_python.py"
try:
    subprocess.run([
        sys.executable, python_script,
        forward_file, reverse_file,
        reads_folder, plots_folder, results_folder
    ], check=True)
    print("\n✅ Python DADA2 pipeline complete!\n")
except subprocess.CalledProcessError as e:
    print(f"\n❌ Python pipeline failed with error code {e.returncode}")
    sys.exit(1)

print("=" * 60)
print("✅ Both DADA2 pipelines complete!")
print("=" * 60)
print(f"📂 Filtered reads saved in: {reads_folder}")
print(f"📊 Quality plots saved in: {plots_folder}")
print(f"📄 Results saved in: {results_folder}")
print(f"\n💡 Run comparison.py to compare R vs Python outputs")
