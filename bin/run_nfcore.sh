NXF_SINGULARITY_CACHEDIR=nxf_singularity_cache \
  nextflow run nf-core/rnaseq \
    -r 3.18.0 \
    -params-file /params.yaml \
    --outdir nf-core-results \
    -ansi-log false \
    -resume
