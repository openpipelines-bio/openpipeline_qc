
import mudata as mu
import pandas as pd

def inspect_h5mu(filepath, name):
    print(f"--- Inspecting {name} ({filepath}) ---")
    try:
        mdata = mu.read_h5mu(filepath)
        rna = mdata.mod['rna']
        print(f"Obs columns ({len(rna.obs.columns)}):")
        print(rna.obs.columns.tolist())
        
        print(f"\nObs head:")
        print(rna.obs.head())
        
        print(f"\nObsm keys: {list(rna.obsm.keys())}")
        
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
    print("\n")

inspect_h5mu("resources_test/spatial_qc_sample_data/visium_tiny.qc.h5mu", "Visium")
