# Workflow template

This repository contains a [nextflow](https://www.nextflow.io/) workflow
template that can be used as the basis for creating new workflows.

> This workflow is not intended to be used by end users.


## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker of conda is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-template --help
```

**Download demo data**

```
# download clair demo data
OUTPUT=output
nextflow run epi2me-labs/wf-clair \
   -w {OUTPUT}/workspace/
   --download 
```

```

# run the pipeline with the test data
OUTPUT=output
nextflow run epi2me-labs/wf-clair \
    -w ${OUTPUT}/workspace \
    -profile standard \
    --bam ${OUTPUT}/demo_data/HG003_chr20_demo.bam \
    --bed ${OUTPUT}/demo_data/quick_demo.bed \
    --ref ${OUTPUT}/demo_data/GRCh38_no_alt_chr20.fa \
    --out_dir ${OUTPUT}
```

The output of the pipeline will be found in `./output` for the above
example. This directory contains the nextflow working directories alongside
the 

TODO: More on outputs

**Workflow outputs**

The primary outputs of the workflow include:

* a simple text file providing a summary of sequencing reads,
* an HTML report document detailing the primary findings of the workflow.


## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)