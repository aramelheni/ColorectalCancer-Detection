from Bio import SeqIO

# Paths to files (Update these to your actual file locations)
dada2_output_path = r"C:\Users\arame\Downloads\SRR29923448_F_filt.fastq"
python_output_path = r"C:\Users\arame\OneDrive\Desktop\Kraya\ISS\CRCProject\filtered\SRR29923448_1_F_filt.fastq"

# Function to read FASTQ into list of (id, sequence, quality)
def read_fastq(file_path):
    records = []
    with open(file_path, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            records.append((record.id, str(record.seq), record.letter_annotations["phred_quality"]))
    return records

# Load both files
dada2_records = read_fastq(dada2_output_path)
python_records = read_fastq(python_output_path)

# Compare number of sequences
if len(dada2_records) != len(python_records):
    print(f" Different number of sequences! DADA2: {len(dada2_records)}, Python: {len(python_records)}")
else:
    print(f" Same number of sequences: {len(dada2_records)}")

# Compare sequences and quality scores
all_match = True
for (dada2_id, dada2_seq, dada2_qual), (py_id, py_seq, py_qual) in zip(dada2_records, python_records):
    if dada2_seq != py_seq:
        print(f"Sequence mismatch in {dada2_id}")
        all_match = False
        break
    if dada2_qual != py_qual:
        print(f" Quality score mismatch in {dada2_id}")
        all_match = False
        break

if all_match:
    print("✅ All sequences and quality scores match exactly!")

else:
    print("⚠️ There were differences between the DADA2 and Python outputs.")
