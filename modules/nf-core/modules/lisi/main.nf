
def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


process LISI{
    // Takes a list of reduced_dims and calculates lisi
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    // label 'process_high'
    label 'process_medium'
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("normalized_pca.h5ad")) {
                        null
                    } else if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
                        null
                    } else {
                        filename.replaceAll("-", "")
                    }
                },
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__metadata)
        val(variables)
        file(file__reduced_dims)
        //tuple(val(label__reduced_dims), file(file__reduced_dims))

    output:
        val(outdir, emit: outdir)
        path(file__metadata, emit: metadata)
        path("${outfile}-lisi.tsv.gz", emit: clusters)
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true

    script:
        outdir = "${outdir_prev}"

        // from the file__metadata job.
        outfile = "${file__metadata}".minus(".tsv.gz")
            .split("-").drop(1).join("-")
        file__reduced_dims = file__reduced_dims.join("::")
        label__reduced_dims = file__reduced_dims
            .replaceAll("reduced_dims-", "")
            .replaceAll(".tsv.gz", "")

        """
        rm -fr plots
        0047-lisi.py \
            --reduced_dims_tsv ${file__reduced_dims} \
            --reduced_dims_tsv_labels ${label__reduced_dims} \
            --metadata_tsv ${file__metadata} \
            --metadata_columns ${variables} \
            --perplexity 30 \
            --output_file ${outfile}-lisi
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """

}
