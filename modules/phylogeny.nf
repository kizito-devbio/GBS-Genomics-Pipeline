process PHYLOGENY {
    tag "phylogeny"
    label 'high_compute'

    conda 'bioconda::iqtree=2.2.0'
    container 'kizitodevbio/strepto-pipeline:latest'

    publishDir "${params.outdir}/phylogeny", mode: 'copy'
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: 'step_log.txt'

    input:
    path aligned

    output:
    path "tree.nwk",      emit: tree
    path "step_log.txt",  emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    TIMESTAMP=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "[\$TIMESTAMP] Step: Phylogenetic Tree Construction using IQ-TREE started." > step_log.txt

    # Run IQ-TREE with Model Finder Plus and bootstrap
    iqtree2 -s "$aligned" -m MFP -bb 1000 -nt AUTO --redo --prefix phylogeny_out >> step_log.txt 2>&1

    # Rename output tree to standard name
    if [ -f "phylogeny_out.treefile" ]; then
        cp "phylogeny_out.treefile" "tree.nwk"
    else
        echo "[\$TIMESTAMP] ERROR: IQ-TREE did not produce a treefile." >> step_log.txt
        touch tree.nwk
        exit 1
    fi

    # Count number of tips in the tree using your specific logic
    num_tips=\$(grep -o "[^,()]*:[0-9]" tree.nwk | wc -l || echo "0")
    echo "[\$TIMESTAMP] Phylogenetic tree constructed with \$num_tips taxa." >> step_log.txt
    echo "[\$TIMESTAMP] PHYLOGENY step completed successfully." >> step_log.txt
    """
}
