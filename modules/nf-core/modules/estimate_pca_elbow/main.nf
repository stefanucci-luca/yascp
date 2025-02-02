
def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}

process ESTIMATE_PCA_ELBOW {
    // Takes annData object, estiamtes the elbow in PC var explained.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("-", "")},
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        val(add_n_to_estimate)

    output:
        val(outdir, emit: outdir)
        path("${outfile}.tsv", emit: pca_elbow_estimate)
        env(AUTO_ELBOW, emit: auto_elbow)
        path("plots/*.png")
        path("plots/*.pdf") optional true

    script:
        
        outdir = "${outdir_prev}"
        log.info("""outdir = ${outdir}""")
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad")
            .split("-").drop(1).join("-")
        outfile = "${outfile}-knee"
        """
            rm -fr plots
            0030-estimate_pca_elbow.py \
                --h5_anndata ${file__anndata} \
                --add_n_pcs_to_elbow ${add_n_to_estimate} \
                --output_file ${outfile}
            mkdir plots
            mv *pdf plots/ 2>/dev/null || true
            mv *png plots/ 2>/dev/null || true
            AUTO_ELBOW=\$(cat ${outfile}-auto_elbow_estimate.tsv)
        """
}
