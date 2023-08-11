process REMOVE_DUPLICATED_DONORS_FROM_GT{
    label 'process_tiny'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }

  input:
    tuple val(samplename),path(subset_genotype),path(subset_genotype_csi)
    path(bridge_file)
    path(input_data_table)

  output:
    // tuple val(pool_id), path("${vireo_fixed_vcf}"), path("${vireo_fixed_vcf}.tbi"), emit: gt_pool
     tuple val(samplename),path("dubs_removed__${subset_genotype}"),path("dubs_removed__${subset_genotype}.csi"),emit:merged_expected_genotypes
    script:
  """
    echo ${bridge_file} > output_vireo.csv
    echo ${samplename}
    echo ${subset_genotype}

    bcftools query -l ${subset_genotype} > donors.tsv
    remove_duplicated_genotypes.py -d donors.tsv -b ${bridge_file} -i ${input_data_table} -n ${samplename}
    bcftools view -S t.tsv ${subset_genotype} -Oz -o dubs_removed__${subset_genotype}
    bcftools index dubs_removed__${subset_genotype}
    rm t.tsv
  """
}

process CAPTURE_VIREO{
  label 'process_tiny'
  input:
    path(vireo_location)
   
  output:
    // tuple val(pool_id), path("${vireo_fixed_vcf}"), path("${vireo_fixed_vcf}.tbi"), emit: gt_pool
    path("output_vireo.csv"),emit:vireo_loc
  script:
  """
    for OUTPUT in \$(ls -d ${vireo_location}/*/)
    do
    samplename1=\$(echo \$OUTPUT | sed 's/${vireo_location}\\///g')
    samplename1=\${samplename1:0:-1}
    echo "\$samplename1 \$PWD/${vireo_location}/\$samplename1/\${samplename1}_headfix_vireo.vcf.gz \$PWD/${vireo_location}/\$samplename1/\${samplename1}_headfix_vireo.vcf.gz.tbi" >> output_vireo.csv
    done
  """    
}



process VIREO_SUBSAMPLING {
    // This module is used to make sure that no cells that there are no cells assigned to the wrong donor.
    // We subsample the cellsnp files to the 80% of random SNPs and run vireo with this.

    tag "${samplename}"
    label 'medium_cpus'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }

    input:
      tuple val(samplename), path(cell_data), val(n_pooled), path(donors_gt_vcf), path(donor_gt_csi), val(itteration)

    output:
      tuple val(samplename), path("vireo_${samplename}___${itteration}"), emit: output_dir
      tuple val(samplename), path("vireo_${samplename}___${itteration}/donor_ids.tsv"), emit: sample_donor_ids
      tuple val(samplename), path("vireo_${samplename}___${itteration}/GT_donors.vireo.vcf.gz"), path(vcf_file),path(donor_gt_csi), emit: sample_donor_vcf
      tuple val(samplename), path("vireo_${samplename}___${itteration}/GT_donors.vireo.vcf.gz"), emit: infered_vcf
      path("vireo_${samplename}___${itteration}/${samplename}.sample_summary.txt"), emit: sample_summary_tsv
      path("vireo_${samplename}___${itteration}/${samplename}__exp.sample_summary.txt"), emit: sample__exp_summary_tsv
      tuple  val(samplename), path("vireo_${samplename}___${itteration}/GT_donors.vireo.vcf.gz"), path("vireo_${samplename}___${itteration}/${samplename}.sample_summary.txt"),path("vireo_${samplename}___${itteration}/${samplename}__exp.sample_summary.txt"),path("vireo_${samplename}___${itteration}/donor_ids.tsv"),path(vcf_file),path(donor_gt_csi), emit: all_required_data
      tuple val(samplename), path("sub_${samplename}_Expected.vcf.gz"), emit: exp_sub_gt optional true
    script:
      vcf_file = ""
      if (params.genotype_input.vireo_with_gt){
        // vcf = " -d ${donors_gt_vcf} --forceLearnGT"
        vcf_file = donors_gt_vcf
        vcf = " -d sub_${samplename}_Expected.vcf.gz --forceLearnGT"
        subset = "bcftools view ${donors_gt_vcf} -R ${cell_data}/cellSNP.cells.vcf.gz -Oz -o sub_${samplename}_Expected.vcf.gz"
        com2 = "cd vireo_${samplename} && ln -s ../${donors_gt_vcf} GT_donors.vireo.vcf.gz"
        com2 = ""
        subset_processing="""
            bcftools view ${donors_gt_vcf} -R ${cell_data}/cellSNP.base.vcf.gz  -Oz -o pre_Overlapping.vcf.gz
            bcftools view -G ${cell_data}/cellSNP.cells.vcf.gz -Oz -o pre_cellSNP.cells.vcf.gz 
            bcftools sort pre_cellSNP.cells.vcf.gz -Oz -o pre_cellSNP3.cells.vcf.gz
            bcftools index pre_cellSNP3.cells.vcf.gz
            bcftools view pre_cellSNP3.cells.vcf.gz -R pre_Overlapping.vcf.gz -Oz -o pre_cellSNP4.cells.vcf.gz
            random_variants.py --random_state ${itteration} --vcf pre_cellSNP4.cells.vcf.gz --rate 60 --order ${cell_data}/cellSNP.base.vcf.gz
            zcat ${cell_data}/cellSNP.base.vcf.gz | head -n2 > subset_${params.vireo.rate}/cellSNP.base.vcf && echo 'next'
            cat random_variants.tsv >> subset_${params.vireo.rate}/cellSNP.base.vcf
            gzip subset_${params.vireo.rate}/cellSNP.base.vcf

            cp ${cell_data}/cellSNP.samples.tsv subset_${params.vireo.rate}/
            # Update the coordinates matrix
            cellsnp_update.R ${cell_data} ./subset_${params.vireo.rate} ./subset_${params.vireo.rate}/cellSNP.base.vcf.gz

        """

      }else{
         vcf = ""
         vcf_file = donors_gt_vcf
         com2 = ""
         subset=""
         reference_expansion_with_piled_up_positions = ""

          subset_processing="""
              # Produce 80% of SNPs panel 
            random_variants.py --random_state ${itteration} --vcf ${cell_data}/cellSNP.base.vcf.gz --rate ${params.vireo.rate}

            # Subset the VCF file
            
            zcat ${cell_data}/cellSNP.base.vcf.gz | head -n2 > subset_${params.vireo.rate}/cellSNP.base.vcf && echo 'next'
            cat random_variants.tsv >> subset_${params.vireo.rate}/cellSNP.base.vcf
            gzip subset_${params.vireo.rate}/cellSNP.base.vcf
            cp ${cell_data}/cellSNP.samples.tsv subset_${params.vireo.rate}/
            # Update the coordinates matrix
            cellsnp_update.R ${cell_data} ./subset_${params.vireo.rate} ./subset_${params.vireo.rate}/cellSNP.base.vcf.gz
              
          """

      }

    """

        mkdir subset_${params.vireo.rate}
        ${subset_processing}
        ${subset}
        umask 2 # make files group_writable
        vireo -c ./subset_${params.vireo.rate} -N $n_pooled -o vireo_${samplename}___${itteration} ${vcf} -t GT --randSeed 1 -p $task.cpus --nInit 200
        # add samplename to summary.tsv,
        # to then have Nextflow concat summary.tsv of all samples into a single file:
        gzip vireo_${samplename}___${itteration}/GT_donors.vireo.vcf || echo 'vireo_${samplename}___${itteration}/GT_donors.vireo.vcf already gzip'
        cat vireo_${samplename}___${itteration}/summary.tsv | \\
          tail -n +2 | \\
          sed s\"/^/${samplename}\\t/\"g > vireo_${samplename}___${itteration}/${samplename}.sample_summary.txt

        cat vireo_${samplename}___${itteration}/summary.tsv | \\
          tail -n +2 | \\
          sed s\"/^/${samplename}__/\"g > vireo_${samplename}___${itteration}/${samplename}__exp.sample_summary.txt
        ${com2}

    """
}


process VIREO {
    tag "${samplename}"
    label 'medium_cpus'
    publishDir "${params.outdir}/deconvolution/vireo/${samplename}/",  mode: "${params.vireo.copy_mode}", overwrite: true,
	  saveAs: {filename -> filename.replaceFirst("vireo_${samplename}/","") }



    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc_with_GT/scrna_deconvolution_v3.img"
        //// container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }

     when:
      params.vireo.run

    input:
      tuple val(samplename), path(cell_data), val(n_pooled), path(donors_gt_vcf), path(donor_gt_csi)

    output:
      tuple val(samplename), path("vireo_${samplename}"), emit: output_dir
      tuple val(samplename), path("vireo_${samplename}/donor_ids.tsv"), emit: sample_donor_ids
      tuple val(samplename), path("vireo_${samplename}/GT_donors.vireo.vcf.gz"), path(vcf_file),path(donor_gt_csi), emit: sample_donor_vcf
      tuple val(samplename), path("vireo_${samplename}/GT_donors.vireo.vcf.gz"), emit: infered_vcf
      path("vireo_${samplename}/${samplename}.sample_summary.txt"), emit: sample_summary_tsv
      path("vireo_${samplename}/${samplename}__exp.sample_summary.txt"), emit: sample__exp_summary_tsv
      tuple  val(samplename), path("vireo_${samplename}/GT_donors.vireo.vcf.gz"), path("vireo_${samplename}/${samplename}.sample_summary.txt"),path("vireo_${samplename}/${samplename}__exp.sample_summary.txt"),path("vireo_${samplename}/donor_ids.tsv"),path(vcf_file),path(donor_gt_csi), emit: all_required_data
      tuple val(samplename), path("sub_${samplename}_Expected.vcf.gz"), emit: exp_sub_gt optional true
    script:
      vcf_file = ""
      if (params.genotype_input.vireo_with_gt){
        vcf = " -d sub_${samplename}_Expected.vcf.gz --forceLearnGT"
        subset = "bcftools view ${donors_gt_vcf} -R ${cell_data}/cellSNP.cells.vcf.gz -Oz -o sub_${samplename}_Expected.vcf.gz"
        vcf_file = donors_gt_vcf
        com2 = "cd vireo_${samplename} && ln -s ../${donors_gt_vcf} GT_donors.vireo.vcf.gz"
        com2 = ""
        // We need to make sure that the genotyes are only informative and not limmiting - for this reason we add in a subset all the variat sites that vireo currently doesnt contain.
        reference_expansion_with_piled_up_positions = "bcftools view ${cell_data}/cellSNP.cells.vcf.gz -G -Oz -o cellsp_piled_up_sites.vcf.gz && bcftools sort cellsp_piled_up_sites.vcf.gz -Oz -o cellsp_piled_up_sites_srt.vcf.gz && bcftools index cellsp_piled_up_sites_srt.vcf.gz && bcftools index sub_${samplename}_Expected.vcf.gz && bcftools merge cellsp_piled_up_sites_srt.vcf.gz sub_${samplename}_Expected.vcf.gz -Oz -o sub_Expected.vcf.gz"
      }else{
         vcf = ""
         vcf_file = donors_gt_vcf
         com2 = ""
         subset=""
         reference_expansion_with_piled_up_positions = ""
      }

    """

      umask 2 # make files group_writable

      ${subset}
      vireo -c $cell_data -N $n_pooled -o vireo_${samplename} ${vcf} -t GT --randSeed 1 -p $task.cpus --nInit 200
      # add samplename to summary.tsv,
      # to then have Nextflow concat summary.tsv of all samples into a single file:
      gzip vireo_${samplename}/GT_donors.vireo.vcf || echo 'vireo_${samplename}/GT_donors.vireo.vcf already gzip'
      cat vireo_${samplename}/summary.tsv | \\
        tail -n +2 | \\
        sed s\"/^/${samplename}\\t/\"g > vireo_${samplename}/${samplename}.sample_summary.txt

      cat vireo_${samplename}/summary.tsv | \\
        tail -n +2 | \\
        sed s\"/^/${samplename}__/\"g > vireo_${samplename}/${samplename}__exp.sample_summary.txt
      ${com2}
    """
}

process VIREO_SUBSAMPLING_PROCESSING{
    tag "${samplename}"
    label 'medium_cpus'
    publishDir  path: "${params.outdir}/concordances/${samplename}",
                mode: "${params.copy_mode}",
                overwrite: "true"


    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }
    input:
      tuple val(samplename), path(vireo_subsampling_folders)

    output:
      tuple val(samplename), path("${samplename}_subsampling_donor_swap_quantification.tsv"), emit: subsampling_donor_swap

    script:   
    """
      echo ${samplename}
      echo ${vireo_subsampling_folders}
      ln -s $projectDir/bin/fix_vireo_header.sh ./fix_vireo_header.sh
      gt_check_and_report_cell_swaps.py
      ln -s subsampling_donor_swap_quantification.tsv ${samplename}_subsampling_donor_swap_quantification.tsv
    """


}

