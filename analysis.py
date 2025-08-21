import pandas as pd
import numpy as np
import os

# Base path setup
base_path = os.path.dirname(__file__)

# Paths to all files (now .csv)
seqtab_path = os.path.join(base_path, 'filtered', 'results', 'seqtab_nochim.csv')
metadata_path = os.path.join(base_path, 'filtered', 'metadata.csv')
taxa_path = os.path.join(base_path, 'filtered', 'results', 'taxa.csv')

# ✅ Load sequence table (CSV)
df_seqtab = pd.read_csv(seqtab_path, index_col=0)
print(f"✅ Sequence Table Loaded: {df_seqtab.shape}")

# ✅ Load metadata (CSV)
metadata = pd.read_csv(metadata_path)
print(f"✅ Metadata Loaded: {metadata.shape}")

# ✅ Load taxonomy table (CSV)
df_taxa = pd.read_csv(taxa_path, index_col=0)
print(f"✅ Taxonomy Table Loaded: {df_taxa.shape}")

# 🧐 Preview loaded data (for debugging)
print("\n🔎 Preview of Sequence Table:")
print(df_seqtab.head())

print("\n🔎 Preview of Metadata:")
print(metadata.head())

print("\n🔎 Preview of Taxonomy Table:")
print(df_taxa.head())

# ✅ Analysis-ready
print("\n✅ All files loaded successfully. Ready for analysis!")
