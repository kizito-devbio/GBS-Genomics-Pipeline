process INTEGRATION_VISUALIZATION {
    tag "visualization"
    label 'low_compute'

    conda 'bioconda::ete3=3.1.3 bioconda::pandas'
    container 'kizitodevbio/strepto-pipeline:latest'

    publishDir "${params.outdir}/viz", mode: 'copy'

    input:
    path tree
    path "inputs/*" // This will receive a flat list of all .tsv files

    output:
    path "metadata.tsv",       emit: metadata
    path "tree_with_meta.nwk", emit: tree_out
    path "tree_viz.png",       emit: png
    path "viz.log",            emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="viz.log"
    TIMESTAMP=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "[\$TIMESTAMP] Integration started." > "\$LOG"

    # Header for the metadata table
    echo -e "SampleID\\tAMR_Profile\\tVF_Profile" > metadata.tsv

    # Loop through the files we collected
    # We use the AMR file as the anchor
    for amr_file in inputs/*_amr.tsv; do
        if [ -f "\$amr_file" ]; then
            # Extract GBS_A909 from inputs/GBS_A909_amr.tsv
            filename=\$(basename "\$amr_file")
            sample=\${filename%_amr.tsv}
            
            # Look for the corresponding VF file in the same inputs folder
            vf_file="inputs/\${sample}_vf.tsv"

            echo "Processing sample: \$sample" >> "\$LOG"

            # Extract AMR genes
            AMR_LIST=\$(awk -F'\\t' 'NR>1 {print \$6}' "\$amr_file" | sort | uniq | paste -sd "," -)
            AMR_LIST=\${AMR_LIST:-None}

            # Extract VF genes if file exists
            if [ -f "\$vf_file" ]; then
                VF_LIST=\$(awk -F'\\t' 'NR>1 {print \$6}' "\$vf_file" | sort | uniq | paste -sd "," -)
                VF_LIST=\${VF_LIST:-None}
            else
                VF_LIST="None"
                echo "Warning: No VF file found for \$sample" >> "\$LOG"
            fi

            echo -e "\${sample}\\t\${AMR_LIST}\\t\${VF_LIST}" >> metadata.tsv
        fi
    done

    # Prepare tree
    cp "$tree" tree_with_meta.nwk

    # Python Visualization logic
    python3 - <<EOF >> "\$LOG" 2>&1
from ete3 import Tree, TreeStyle, NodeStyle, TextFace
import pandas as pd
import os

t = Tree("tree_with_meta.nwk", format=1)
df = pd.read_csv("metadata.tsv", sep="\\t")

# Create lookup dictionaries
amr_map = df.set_index('SampleID')['AMR_Profile'].to_dict()
vf_map = df.set_index('SampleID')['VF_Profile'].to_dict()

for leaf in t.iter_leaves():
    # Clean name in case of .fna or .fasta suffixes in the tree
    clean_name = leaf.name.replace(".fna", "").replace(".fasta", "").replace(".fsa", "")
    
    if clean_name in amr_map:
        nstyle = NodeStyle()
        nstyle['fgcolor'] = 'red' if amr_map[clean_name] != 'None' else 'green'
        nstyle['size'] = 12
        leaf.set_style(nstyle)
        
        # Add labels for Genes
        amr_text = f" [AMR: {amr_map[clean_name]}]"
        vf_text = f" [VF: {vf_map[clean_name]}]"
        
        leaf.add_face(TextFace(amr_text, fsize=10, fgcolor="red"), column=1, position="branch-right")
        leaf.add_face(TextFace(vf_text, fsize=10, fgcolor="blue"), column=2, position="branch-right")

ts = TreeStyle()
ts.show_leaf_name = True
ts.title.add_face(TextFace("Streptococcus Pangenome Analysis", fsize=16), column=0)
t.render("tree_viz.png", tree_style=ts)
EOF

    echo "[\$TIMESTAMP] Visualization finished." >> "\$LOG"
    """
}
