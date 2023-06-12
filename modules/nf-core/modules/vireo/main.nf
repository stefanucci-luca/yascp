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
    label 'process_medium'

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
        vcf = " -d sub_${samplename}_Expected.vcf.gz --forceLearnGT"
        subset = "bcftools view ${donors_gt_vcf} -R ${cell_data}/cellSNP.cells.vcf.gz -Oz -o sub_${samplename}_Expected.vcf.gz"
        vcf_file = donors_gt_vcf
        com2 = "cd vireo_${samplename}___${itteration} && ln -s ../${donors_gt_vcf} GT_donors.vireo.vcf.gz"
        com2 = ""
      }else{
         vcf = ""
         vcf_file = donors_gt_vcf
         com2 = ""
         subset=""
      }

    """

      # Produce 80% of SNPs panel 
      random_variants.py --random_state ${itteration} --vcf ${cell_data}/cellSNP.base.vcf.gz

      # Subset the VCF file
      mkdir subset_80
      zcat ${cell_data}/cellSNP.base.vcf.gz | head -n2 > subset_80/cellSNP.base.vcf && echo 'next'
      cat random_variants.tsv >> subset_80/cellSNP.base.vcf
      gzip subset_80/cellSNP.base.vcf
      cp ${cell_data}/cellSNP.samples.tsv subset_80/
      # Update the coordinates matrix
      cellsnp_update.R ${cell_data} ./subset_80 ./subset_80/cellSNP.base.vcf.gz

      umask 2 # make files group_writable
      ${subset}
      vireo -c ./subset_80 -N $n_pooled -o vireo_${samplename}___${itteration} ${vcf} -t GT --randSeed 1 -p $task.cpus --nInit 200
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


process VIREO_SUBSAMPLING_PROCESSING{
    tag "${samplename}"
    label 'process_medium'
    publishDir "${params.outdir}/deconvolution/vireo/${samplename}/",  mode: "${params.vireo.copy_mode}", overwrite: true


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

process VIREO {
    tag "${samplename}"
    label 'process_high'
    publishDir "${params.outdir}/deconvolution/vireo/${samplename}/",  mode: "${params.vireo.copy_mode}", overwrite: true,
	  saveAs: {filename -> filename.replaceFirst("vireo_${samplename}/","") }



    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
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
      }else{
         vcf = ""
         vcf_file = donors_gt_vcf
         com2 = ""
         subset=""
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
