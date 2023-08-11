// Default parameters for bbknn
params.bbknn = [
    run_process: false,
    batch_variable: [value: ["experiment_id"]]
]

def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


process BBKNN{
    
    // Calulates bbknn neighbors and saves UMAPS of these
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
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
                    } else if(filename.endsWith("${param_details}.tsv.gz")) {
                        null
                    } else {
                        filename.replaceAll("-", "")
                    }
                },
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__metadata)
        path(file__pcs)
        val(params__pcs)
        each n_pcs
        each batch_var

    output:
        val(outdir, emit: outdir)
        path("${outfile}-bbknn.h5ad", emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path("reduced_dims.tsv.gz", emit: reduced_dims)
        // NOTE: passing the param details as an unpublished file is messy,
        // but I could not get collect of ${param_details} and file to work.
        path(
            "reduced_dims-${param_details}.tsv.gz",
            emit: reduced_dims_params
        )

    script:
        param_details = "${params__pcs}-bbknn"
        param_details = "${param_details}.batch=${batch_var}"
        param_details = "${param_details}.n_pcs=${n_pcs}"
        outdir = "${outdir_prev}/reduced_dims-${param_details}"
        // from the file__anndata job.
        outfile = "outfile_adata"
        """
        
        0045-bbknn.py \
            --h5_anndata ${file__anndata} \
            --batch_key ${batch_var} \
            --n_pcs ${n_pcs} \
            --output_file ${outfile}-bbknn
        cp ${outfile}-bbknn-reduced_dims.tsv.gz \
            reduced_dims-${param_details}.tsv.gz
        mv ${outfile}-bbknn-reduced_dims.tsv.gz \
            reduced_dims.tsv.gz
        """

}
