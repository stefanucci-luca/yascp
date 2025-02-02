process KERAS_CELLTYPE {

    tag { "${experiment_id}" }
    publishDir  path: "${params.outdir}/celltype/keras_celltype/${experiment_id}/",
        mode: "${params.copy_mode}",
        overwrite: "true"

    label 'process_high_memory'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"

    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }
    input:
        tuple(
            val(experiment_id),
            path(keras_input_h5ad) // from cellbender into scrublet
        )
    
    output:
        tuple(
            val(experiment_id),
            path("*.png"),
            path("*_gtr_*.h5ad"),// do not capture input Pilot_study_of_dissociation_methods_for_human_gut_tissues7980357-scrublet.h5ad
            path("*predictions.h5ad"), 
            emit: out1)
        tuple(
            val(experiment_id),
            path("*predictions.h5ad"), 
            emit: to_merge)
        path('*_celltypes.tsv', emit:predicted_celltype_labels)
    
    script:
        """
            0057-predict_clusters_keras_model-anndata.py \\
            --h5_anndata \"${keras_input_h5ad}\" \\
            --h5_layer \"${params.celltype_prediction.keras.h5_layer}\" \\
            --keras_model \"${params.celltype_prediction.keras.keras_model}\" \\
            --keras_weights_df \"${params.celltype_prediction.keras.keras_weights_df}\" \\
            --keras_model_cluster_labels \"${params.celltype_prediction.keras.keras_model_cluster_labels}\" \\
            --filter_top_cell_probabilities \"${params.celltype_prediction.keras.filter_top_cell_probabilities}\" \\
            --output_file \"${experiment_id}___cellbender_fpr${params.cellbender_resolution_to_use}-scrublet-ti_freeze003_prediction\" ${params.celltype_prediction.keras.save_all_probabilities} 
        """
}
