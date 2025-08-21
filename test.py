import os
import subprocess

# Get current project folder
project_folder = os.path.dirname(os.path.abspath(__file__))

# Create filtered folder and subfolders
filtered_folder = os.path.join(project_folder, "filtered")
reads_folder = os.path.join(filtered_folder, "reads")
plots_folder = os.path.join(filtered_folder, "plots")
results_folder = os.path.join(filtered_folder, "results")

for folder in [filtered_folder, reads_folder, plots_folder, results_folder]:
    os.makedirs(folder, exist_ok=True)

# Input files
forward_file = r"C:\Users\arame\OneDrive\Desktop\SRA files\SRR29923448_1.fastq"
reverse_file = r"C:\Users\arame\OneDrive\Desktop\SRA files\SRR29923448_2.fastq"

# Path to R script
rscript_path = r"C:\Program Files\R\R-4.4.2\bin\Rscript.exe"
r_script = os.path.join(project_folder, "dataAssessing.R")
silva_folder = os.path.join(project_folder, "silva")

# Run R script with all required paths
subprocess.run([
    rscript_path, r_script,
    forward_file, reverse_file,
    reads_folder, plots_folder, results_folder, silva_folder
], check=True)

print("✅ Full DADA2 processing pipeline complete!")
print(f"📂 Filtered reads saved in: {reads_folder}")
print(f"📊 Quality plots saved in: {plots_folder}")
print(f"📄 Results saved in: {results_folder}")
