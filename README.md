
<p align="center">
  <img src="https://github.com/wtsi-hgi/yascp/blob/main/assets/images/YASCP_Logo.png" width="80%"/>
</p>

[![GitHub Actions CI Status](https://github.com/nf-core/yascp/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/yascp/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/yascp/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/yascp/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/yascp/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23yascp-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/yascp)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

- ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) `Pipeline is currently under develpment and therefore there is no stable release yet.`


## Introduction

**nf-core/yascp** is a bioinformatics best-practice analysis pipeline for deconvolution, qc, clustering of a single cell datasets.
This is a large scale single-cell pipeline developed initially for processing Cardinal project samples, however it is applicable to any other scRNA analysis. The pipeline has been inspired by deconvolution (https://github.com/wtsi-hgi/nf_scrna_deconvolution.git ), cellbender (https://github.com/wtsi-hgi/nf_cellbender ) and qc (https://github.com/wtsi-hgi/nf_qc_cluster/tree/main ) pipelines. Input requires a tsv seperated file with paths to the Cellranger 6.11 outputs (however we will shortly add a Cellranger module to make this pipeline more transferable) and if running in an genotype  additional input is required to be provided in an input.nf file pointing to the vcf location. This pipeline is designed to be used for multiple large scale single cell experiments.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->
<!-- On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/yascp/results). -->

## Pipeline summary

1. Cellbender
2. CellSNP
3. Vireo
4. Souporcell
5. Celltypist
6. Azimuth
7. BBKNN
8. Harmony
9. Scrublet
10. Sccaf
11. Lisi
12. Isolation Forest
13. Hard filters
14. Genotype deconvolution and GT match against multiple panels.
15. Citeseq DSB normalisations
16. Cell genotype concordance Calculations

<p align="center">
  <img src="https://github.com/wtsi-hgi/yascp/blob/main/assets/images/yascp_workflow.png" width="100%"/>
</p>

## Documentation: Prepearing your own data and interpreting the results

To understand how to prepeare your own data and how to interpret the results please refear to documents [HERE](https://github.com/wtsi-hgi/yascp/tree/main/docs)

## Quick Start
Easyest to do is using a conda enviroment.

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`)
    ```console
    conda install nextflow=21.04.0
    ```

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download/clone the pipeline and test it on a minimal dataset with a single command:

    !NOTE: you need to define your institution specific queues in the conf/base.conf or provide aditional config file with -c flag in folowing comand such as: -c /path/to/yascp/conf/extra_confs/sanger/base.conf



    ```console
    nextflow run /path/to/colned/yascp -profile test_full,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

    !ALSO: by default this test dataset will run cellbender with a cpus - only performing 10 epochs. Cellbender is built for a gpu queue, so for the actual runs the deafault is to utilise_gpu = true
    You need to make sure that the gpu queue according to your institution is defined in confs: withLabel: gpu {} as for sanger [`config file`](https://github.com/wtsi-hgi/yascp/blob/0fce7bd8ce4ca734e34b28443fe89630e295b1eb/conf/extra_confs/sanger/base.conf#L180-L194)

    > * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
    > * If you are using `singularity` then the pipeline will auto-detect this and attempt to download the Singularity images directly as opposed to performing a conversion from Docker images. If you are persistently observing issues downloading Singularity images directly due to timeout or network issues then please use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. Alternatively, it is highly recommended to use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to pre-download all of the required containers before running the pipeline and to set the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options to be able to store and re-use the images from a central location for future pipeline runs.



## Credits

Yascp was originally written by Matiss Ozols as part of the Cardinal project with contributions from Leland Taylor, Guillaume Noell, Hannes Ponstingl, Vivek Iyer,  Henry Taylor, Tobi Alegbe.

We thank the following people for their extensive assistance in the development of this pipeline:
Monika Krzak
<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#yascp` channel](https://nfcore.slack.com/channels/yascp) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

Currently pipeline has not been published but we would appreciate if you coul please acknowlage the use of this pipeline in your work.

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/yascp for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->
<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

## T regopathies branch

Analysis pipeline tailored for the analysis of the scRNA-Seq for the T-regopathies project
