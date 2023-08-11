process AZIMUTH{
    tag "${samplename}"    
    label 'process_medium'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_azimuth_d54db9b-2021-12-13-8dd0b7fce918.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/seurat_azimuth_pbmc_1.0.img"
        
    } else {
        container "wtsihgi/nf_scrna_qc_azimuth:d54db9b"
    }

    publishDir  path: "${params.outdir}/celltype/azimuth/",
            saveAs: {filename -> "${outfil_prfx}_" + filename},
            mode: "${params.copy_mode}",
            overwrite: "true"
    stageInMode 'copy'
    // stageInMode 'copy' is needed because SeuratDisk:::Convert()
    // generates the output file apparently from the absolute path of input file.
    // Symbolic links have the output file written to the link target directory
    // where it cannot be found by the azimuth.R script.

    input:
        val outdir_prev
        path file_h5ad_batch

    output:
        path(celltype_table, emit:predicted_celltypes)
        path(celltype_table, emit:predicted_celltype_labels)
        path "ncells_by_type_barplot.pdf"
        path "query_umap.pdf"
        path "prediction_score_umap.pdf"
        path "prediction_score_vln.pdf"
        path "mapping_score_umap.pdf"
        path "mapping_score_vln.pdf"

    script:
    
    outdir = "${outdir_prev}/azimuth"
    // output file prefix: strip random hex number form beginning of file name
    outfil_prfx = "${file_h5ad_batch}".minus(".h5ad")
    //outfil_prfx = "${file_h5ad_batch}".minus(".h5ad")
    celltype_table = "${outfil_prfx}_predicted_celltype_l2.tsv.gz"
    """
        azimuth.R ./${file_h5ad_batch}
        gzip -c predicted_celltype_l2.tsv > ${celltype_table}
    """
}

process REMAP_AZIMUTH{
    // This process remaps Azimuth L2 to L1 and L0
    tag "${samplename}"    
    label 'process_low'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/seurat_azimuth_pbmc_1.0.img"
    } else {
        container "wtsihgi/nf_scrna_qc_azimuth:d54db9b"
    }

    publishDir  path: "${params.outdir}/celltype/azimuth/",
            mode: "${params.copy_mode}",
            overwrite: "true"
    stageInMode 'copy'  

    input:
        path(azimuth_file)
        path(mapping_file)

    output:
        path(celltype_table, emit:predicted_celltype_labels)

    script:
        celltype_table = "remapped__${azimuth_file}"
        """
            remap_azimuth_l2.py -of ${celltype_table} -m ${mapping_file} -az ${azimuth_file}
        """

}
