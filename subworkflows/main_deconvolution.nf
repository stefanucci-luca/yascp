nextflow.enable.dsl=2

// main deconvolution modules, common to all input modes:

include { CELLSNP;capture_cellsnp_files;DYNAMIC_DONOR_EXCLUSIVE_SNP_SELECTION; ASSESS_CALL_RATE } from "$projectDir/modules/nf-core/modules/cellsnp/main"
include { SUBSET_GENOTYPE } from "$projectDir/modules/nf-core/modules/subset_genotype/main"
include { VIREO;REMOVE_DUPLICATED_DONORS_FROM_GT;VIREO_SUBSAMPLING;VIREO_SUBSAMPLING_PROCESSING } from "$projectDir/modules/nf-core/modules/vireo/main"
include { GUZIP_VCF } from "$projectDir/modules/nf-core/modules/guzip_vcf/main"
include { SOUPORCELL } from "$projectDir/modules/nf-core/modules/souporcell/main"
include { SPLIT_DONOR_H5AD } from "$projectDir/modules/nf-core/modules/split_donor_h5ad/main"
include { PLOT_DONOR_CELLS } from "$projectDir/modules/nf-core/modules/plot_donor_cells/main"
include {MULTIPLET} from "$projectDir/modules/nf-core/modules/multiplet/main"
include {SOUPORCELL_VS_VIREO} from "$projectDir/modules/nf-core/modules/plot_souporcell_vs_vireo/main"
include { REPLACE_GT_ASSIGNMENTS_WITH_PHENOTYPE; ENHANCE_VIREO_METADATA_WITH_DONOR } from "$projectDir/modules/nf-core/modules/genotypes/main"
include { match_genotypes } from './match_genotypes'
include {ENHANCE_STATS_GT_MATCH } from "$projectDir/modules/nf-core/modules/genotypes/main"
include {SUBSET_WORKF} from "$projectDir/modules/nf-core/modules/subset_genotype/main"
include {REPLACE_GT_DONOR_ID2; REPLACE_GT_DONOR_ID2 as REPLACE_GT_DONOR_ID_SUBS } from "$projectDir/modules/nf-core/modules/genotypes/main"
include { RETRIEVE_RECOURSES } from "$projectDir/subworkflows/local/retrieve_recourses"
include {GT_MATCH_POOL_IBD } from "$projectDir/modules/nf-core/modules/genotypes/main"
include { MATCH_GT_VIREO } from "$projectDir/modules/nf-core/modules/genotypes/main"

include {VIREO_GT_FIX_HEADER; VIREO_GT_FIX_HEADER as VIREO_GT_FIX_HEADER_SUBS; VIREO_ADD_SAMPLE_PREFIX; MERGE_GENOTYPES_IN_ONE_VCF as MERGE_GENOTYPES_IN_ONE_VCF_INFERED; MERGE_GENOTYPES_IN_ONE_VCF as MERGE_GENOTYPES_IN_ONE_VCF_SUBSET} from "$projectDir/modules/nf-core/modules/genotypes/main"
include {collect_file as collect_file1;
        collect_file as collect_file2;
        collect_file as collect_file3;
        collect_file as collect_file4;
        collect_file as collect_file5;
        collect_file as collect_file6;
        collect_file as collect_file7;
        collect_file as collect_file8;
        collect_file as collect_file9;
        collect_file as collect_file10} from "$projectDir/modules/nf-core/modules/collect_file/main"

workflow  main_deconvolution {

    take:
		ch_experiment_bam_bai_barcodes
		ch_experiment_npooled
		ch_experiment_filth5
		ch_experiment_donorsvcf_donorslist
        channel__file_paths_10x
        vcf_input

    main:
		log.info "#### running DECONVOLUTION workflow #####"

        if (params.reference_assembly_fasta_dir=='https://yascp.cog.sanger.ac.uk/public/10x_reference_assembly'){
            RETRIEVE_RECOURSES()  
            genome = RETRIEVE_RECOURSES.out.reference_assembly
        }else{
            genome = "${params.reference_assembly_fasta_dir}"
        }

        // genome.subscribe { println "genome: $it" }

        if (params.genotype_input.run_with_genotype_input) {
            // We have to produce a single vcf file for each individual pool.
            // Therefore we create 2 channels:
            // 1) All the expected vcf ids listed in the donor table
            Channel.fromPath(params.input_data_table,      
            followLinks: true,
            checkIfExists: true
            ).splitCsv(header: true, sep: '\t').map { row -> tuple(row.experiment_id, row.donor_vcf_ids) }
            .set { donors_in_pools }
            
            // 2) All the vcfs provided to us. 
            vcf_input.splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.label, file(row.vcf_file_path), file("${row.vcf_file_path}.csi")) }
            .set { ch_ref_vcf }

            // This will subsequently result in a joint vcf file per cohort per donors listed for each of the pools that can be used in VIREO and/or GT matching algorythm.
            SUBSET_WORKF(ch_ref_vcf,donors_in_pools,'AllExpectedGT',genome)
            merged_expected_genotypes = SUBSET_WORKF.out.merged_expected_genotypes
            merged_expected_genotypes2 = merged_expected_genotypes.combine(Channel.fromPath(params.cellsnp.vcf_candidate_snps))
            // merged_expected_genotypes2.subscribe { println "merged_expected_genotypes2: $it" }
            GT_MATCH_POOL_IBD(SUBSET_WORKF.out.samplename_subsetvcf_ibd,'Withing_expected','Expected')
            DYNAMIC_DONOR_EXCLUSIVE_SNP_SELECTION(merged_expected_genotypes2)
            cellsnp_panels = DYNAMIC_DONOR_EXCLUSIVE_SNP_SELECTION.out.cellsnp_pool_panel
            informative_uninformative_sites = DYNAMIC_DONOR_EXCLUSIVE_SNP_SELECTION.out.informative_uninformative_sites

            // // If we have selected that we want to use all the genotypes as an input in the VCF file we will use the output of the MERGE_GENOTYPES_IN_ONE_VCF_SUBSET
            // if (params.genotype_input.subset_vireo_genotypes){
            //     MERGE_GENOTYPES_IN_ONE_VCF_SUBSET(SUBSET_WORKF.out.study_merged_vcf.collect(),'subset')
            //     merged_expected_genotypes = MERGE_GENOTYPES_IN_ONE_VCF_SUBSET.out.study_merged_vcf
            // }

        }else{
            ch_ref_vcf = Channel.of()
        }


        
        ch_poolid_donor_assignment = Channel.empty()
        ch_experiment_donorsvcf_donorslist.map { experiment, donorsvcf, donorstbi,donorslist -> tuple(experiment, donorslist.replaceAll(/,/, " ").replaceAll(/"/, ""))}.set{donors_in_lane}
        ch_experiment_bam_bai_barcodes.combine(ch_experiment_npooled, by: 0).set{cellsnp_with_npooled}

        if (params.genotype_input.run_with_genotype_input) {
            // Here we provide sites to pile up the snps within the pool. 
            // As a starting point instead of the default cellsnp panel users are encouraged to:
            if (params.provide_within_pool_donor_specific_sites_for_pilup){
                cellsnp_with_npooled2 = cellsnp_with_npooled.combine(cellsnp_panels, by: 0)
            }else{
                cellsnp_with_npooled2 = cellsnp_with_npooled.combine(Channel.fromPath(params.cellsnp.vcf_candidate_snps))
            }
        }else{
            cellsnp_with_npooled2 = cellsnp_with_npooled.combine(Channel.fromPath(params.cellsnp.vcf_candidate_snps))
        }
        // cellsnp_with_npooled2.subscribe { println "cellsnp_with_npooled2: $it" }
        

        log.info('Capturing some of the existing CELLSNP files')
        capture_cellsnp_files(params.existing_cellsnp)
        capture_cellsnp_files.out.cellsnp_loc.splitCsv(header: false, sep: ' ')
            .map{row->tuple(row[0], "${row[1]}")}
            .set{cellsnp_output_dir1}
        cellsnp_output_dir1.join(cellsnp_with_npooled2, remainder: true).set{filter_channel}
        // filter_channel.subscribe { println "cellsnp_output_dir1 filter_channel: $it" }
        filter_channel.filter{ it[1] == null }.map{row -> tuple(row[0], row[2],row[3],row[4],row[5],row[6])}.set{cellsnp_with_npooled}
        
        capture_cellsnp_files.out.cellsnp_loc.splitCsv(header: false, sep: ' ')
            .map{row->tuple(row[0], file("${row[1]}/cellSNP.cells.vcf.gz"))}
            .set{cellsnp_cell_vcfs}           

        CELLSNP(cellsnp_with_npooled)
        if (params.genotype_input.run_with_genotype_input) {
            // Here we assess how many informative sites has been called on. 
            CELLSNP.out.cell_vcfs.combine(DYNAMIC_DONOR_EXCLUSIVE_SNP_SELECTION.out.informative_uninformative_sites, by: 0).set{assess_call_rate_input}
            ASSESS_CALL_RATE(assess_call_rate_input)
            collect_file9(ASSESS_CALL_RATE.out.variants_description.collect(),"all_variants_description.tsv",params.outdir+'/concordances',1,'')
        }
        


        cellsnp_output_dir2 = CELLSNP.out.cellsnp_output_dir
        cellsnp_output_dir=cellsnp_output_dir1.concat(cellsnp_output_dir2)
        cellsnp_cell_vcfs2=cellsnp_cell_vcfs.concat(CELLSNP.out.cell_vcfs)


        // channel__file_paths_10x.subscribe { println "channel__file_paths_10x: $it" }
        // if (params.sample_qc.cell_filters.filter_multiplets.run_process){
        MULTIPLET(
            channel__file_paths_10x,
            params.sample_qc.cell_filters.filter_multiplets.expected_multiplet_rate,
            params.sample_qc.cell_filters.filter_multiplets.n_simulated_multiplet,
            params.sample_qc.cell_filters.filter_multiplets.multiplet_threshold_method,
            params.sample_qc.cell_filters.filter_multiplets.scale_log10
        )
        scrublet_paths = MULTIPLET.out.scrublet_paths
        // }else{
        //     scrublet_paths = Channel.of()
        // }

//         // Here we run Vireo software to perform the donor deconvolution. Note that we have coded the pipeline to be capable in using
//         // the full genotypes as an input and also subset to the individuals provided as an input in the donor_vcf_ids column.
        if (params.genotype_input.vireo_with_gt && params.genotype_input.run_with_genotype_input) {
            log.info "---running Vireo with genotype input----"
            // for each experiment_id to deconvolute, subset donors vcf to its donors and subset genomic regions.
            // Here we subset the genotypes. This happens if the input.nf contains subset_genotypes = true
            log.info "---We are using subset genotypes running Vireo----"
            // We need to make sure that the expected genotypes dont contain repeated genotypes - donors sequenced twice.
            REMOVE_DUPLICATED_DONORS_FROM_GT(merged_expected_genotypes,params.genotype_phenotype_mapping_file,params.input_data_table)
            merged_expected_genotypes = REMOVE_DUPLICATED_DONORS_FROM_GT.out.merged_expected_genotypes
            cellsnp_output_dir.combine(ch_experiment_npooled, by: 0)
                .combine(merged_expected_genotypes, by: 0).set{full_vcf}
        } else {
            // Vireo can also be run without the genotypes, and the performance is equally good then running with, however have to be careful in what sites file is provided.
            // Here we run it without the genotypes and the donors are labeled as donor0, donor1 etc, dependant on the number of donors set in the input file n_pooled column
            log.info "-----running Vireo without genotype input----"
            cellsnp_output_dir.combine(ch_experiment_npooled, by: 0).set{full_vcf}
            full_vcf.map {experiment, cellsnp, npooled -> tuple(experiment, cellsnp, npooled,[],[])}.set{full_vcf}
        }

        // When all the channels are prpeared accordingly we exacute the vireo with the prpeared channel.
        full_vcf.filter { experiment, cellsnp, npooled, t,ti -> npooled != '1' }.set{full_vcf2}
        full_vcf.filter { experiment, cellsnp, npooled, t,ti -> npooled == '1' }.set{not_deconvoluted}

        //Vireo creates per pool vcf vit donor0 etc, or if genotypes are provided genotype IDs will be placed there.
        // full_vcf2.subscribe { println "full_vcf2: $it" }

        Channel.of(1..params.vireo.subsample_times).set{itterations}

        full_vcf2.combine(itterations).set{vireo_extra_repeats}
        // vireo_extra_repeats.subscribe { println "vireo_extra_repeats is: $it" }
        VIREO(full_vcf2)


        
        //But to make it consistent we still randomly assign donor IDs per pool for each of the names.
        vireo_with_gt = Channel.of(params.genotype_input.vireo_with_gt)
        VIREO.out.all_required_data.set{replacement_input}
        replacement_input.combine(vireo_with_gt).set{vir_repl_input}
        REPLACE_GT_DONOR_ID2(vir_repl_input)
        cell_assignments = REPLACE_GT_DONOR_ID2.out.cell_assignments
        VIREO_GT_FIX_HEADER(REPLACE_GT_DONOR_ID2.out.infered_vcf,genome)
        VIREO_ADD_SAMPLE_PREFIX(VIREO_GT_FIX_HEADER.out.infered_vcf)
        
        // ch_vireo_donor_n_cells_tsv = collect_file9(VIREO_SUBSAMPLING_PROCESSING.out.subsampling_donor_swap.collect(),"vireo_donor_n_cells.tsv",params.outdir+'/deconvolution/filepaths',0,header_seed)

        // VIREO_ADD_SAMPLE_PREFIX.out.infered_vcf.collect().subscribe { println "MERGE_GENOTYPES_IN_ONE_VCF_INFERED: $it" }
        // MERGE_GENOTYPES_IN_ONE_VCF_INFERED(VIREO_ADD_SAMPLE_PREFIX.out.infered_vcf.collect(),'infered')
        
        vireo_out_sample_donor_vcf = VIREO_GT_FIX_HEADER.out.infered_vcf
        vireo_out_sample_summary_tsv = REPLACE_GT_DONOR_ID2.out.sample_summary_tsv
        vireo_out_sample__exp_summary_tsv = REPLACE_GT_DONOR_ID2.out.sample__exp_summary_tsv
        vireo_out_sample_donor_ids = REPLACE_GT_DONOR_ID2.out.sample_donor_ids

        if (params.genotype_input.run_with_genotype_input) {
            VIREO_SUBSAMPLING(vireo_extra_repeats)
            VIREO_SUBSAMPLING.out.output_dir.concat(VIREO.out.output_dir).set{tuple_1}
            tuple_1.groupTuple(by:0).set{vspp0}
            VIREO_SUBSAMPLING_PROCESSING(vspp0)
            VIREO_SUBSAMPLING.out.all_required_data.set{replacement_input_sub}
            // replacement_input_sub.combine(vireo_with_gt).set{vir_repl_input}
            // REPLACE_GT_DONOR_ID_SUBS(vir_repl_input)
            // VIREO_GT_FIX_HEADER_SUBS(REPLACE_GT_DONOR_ID_SUBS.out.infered_vcf,genome)
            // VIREO_GT_FIX_HEADER_SUBS.out.gt_pool
            //     .combine(ch_ref_vcf).set { gt_math_pool_against_panel_input_subs }
                
            // MATCH_GT_VIREO(gt_math_pool_against_panel_input_subs)


            VIREO_GT_FIX_HEADER.out.gt_pool
                .combine(ch_ref_vcf)
                .set { gt_math_pool_against_panel_input }

            match_genotypes(vireo_out_sample_donor_vcf,merged_expected_genotypes,VIREO_GT_FIX_HEADER.out.gt_pool,gt_math_pool_against_panel_input,genome,ch_ref_vcf,cellsnp_cell_vcfs2,cell_assignments,VIREO_SUBSAMPLING_PROCESSING.out.subsampling_donor_swap,informative_uninformative_sites)
            gt_matches = match_genotypes.out.donor_match_table.collect()

            ENHANCE_STATS_GT_MATCH(match_genotypes.out.donor_match_table_enhanced)
            collect_file1(ENHANCE_STATS_GT_MATCH.out.assignments.collect(),"assignments_all_pools.tsv",params.outdir+'/deconvolution/vireo_gt_fix',1,'')
            collect_file10(ENHANCE_STATS_GT_MATCH.out.assignments.collect(),"assignments_all_pools.tsv",params.outdir+'/gtmatch',1,'')
            assignments_all_pools = collect_file1.out.output_collection

        }else{
            gt_matches = Channel.from("$projectDir/assets/fake_file.fq")
            assignments_all_pools = Channel.from("$projectDir/assets/fake_file.fq")
        }

        // if (params.souporcell.run){
        //     // YASCP pipeline is also capable in running SOUPORCELL instead of VIREO. If activated SOUPORCELL will be used.
        //     // yascp currently doesnt have an option to take souporcell assignments as downstream instead of vireo but this will be added shortly.
        //     // This runs the Souporcell Preprocessing

        //     if (params.souporcell.use_raw_barcodes) {
        //         // read raw cellranger barcodes per pool for souporcell
        //         log.info """Here use the raw barcodes"""
        //     } else if (! params.souporcell.use_raw_barcodes) {
        //         ch_experiment_bam_bai_barcodes
        //             .combine(ch_experiment_npooled, by: 0)
        //             .set{ch_experiment_bam_bai_barcodes_npooled}
        //     }

        //     // This runs the Souporcell
        //     // Similarly to the VIREO Soupocell can be run with and without genotypes and folowing prpeares the inputs accordingly to each option.
        //     // Soupocell cant digest a gz file gence we extract the data.
        //     if (params.genotype_input.run_with_genotype_input) {

        //         if (params.genotype_input.subset_genotypes){
        //             // this will run the soupocell with the subset genotypes. This happens if the input.nf contains subset_genotypes = true
        //             log.info "---We are using subset genotypes running Suporcell----"
        //             GUZIP_VCF(SUBSET_GENOTYPE.out.samplename_subsetvcf)
        //             // ch_experiment_bam_bai_barcodes_npooled.combine(SUBSET_GENOTYPE.out.samplename_subsetvcf, by: 0).set{full_vcf}
        //             ch_experiment_bam_bai_barcodes_npooled.combine(GUZIP_VCF.out.souporcell_vcf, by: 0).set{full_vcf}

        //         }else{
        //             log.info "---We are using a full genotype input for Suporcell----"
        //             // this however currently doesnt work and the individuals have to be provided.
        //             // here just add the full vcf path to each of the ch_experiment_bam_bai_barcodes_npooled
        //             GUZIP_VCF(tuple('full_vcf', file(params.genotype_input.full_vcf_file)))
        //             GUZIP_VCF.out.souporcell_vcf.map { sample, vcf -> vcf }.set{vcf_file}
        //             ch_experiment_bam_bai_barcodes_npooled.map {
        //                 samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters ->
        //                 tuple(samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters)
        //                 }.set{full_vcf}
        //             full_vcf.combine(vcf_file).set{full_vcf}

        //         }
        //         full_vcf.combine(donors_in_lane, by: 0).set{full_vcf}
        //     }
        //     else{
        //         // Soupocell can also run without the genotypes. This will prpeare channels to run it withoutt.
        //         log.info "-----running Suporcell without genotype input----"
        //         // here make add an empty entry [] to the ch_experiment_bam_bai_barcodes_npooled
        //         ch_experiment_bam_bai_barcodes_npooled.map {
        //             samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters ->
        //             tuple(samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters,[],[])
        //             }.set{full_vcf}
        //     }

        //     // When all the channels are prepeared then we can run the soupocell accordingly.
        //     SOUPORCELL(full_vcf,
        //         Channel.fromPath(params.souporcell.reference_fasta).collect())
        //     // Regardless if the Soupocell is run with or without genotypes we still need to match the donor ids with the cluster ids since this does not happen automatically is Soupocell.
        // }

        // if (params.souporcell.run && params.vireo.run) {
        //     // This is an orringinal comparison of Vireo and soupocell, this needs to be cleaned up.
        //     SOUPORCELL_VS_VIREO(
        //         vireo_out_sample_donor_ids
        //         .combine(SOUPORCELL.out.souporcell_output_files.map {a,b,c,d -> tuple(a,b)},
        //             by: 0))
        // }


        ch_experiment_bam_bai_barcodes
            .map { samplename, bam_file, bai_file, barcodes_tsv_gz -> tuple(samplename, file(bam_file)) }
            .combine(vireo_out_sample_donor_ids, by: 0 )
            .set { ch_experiment_bam_vireo_donor_ids }

        // If sample is not deconvoluted we will use scrublet to detect the doublets and remove them.
        not_deconvoluted.map{ experiment, donorsvcf, npooled,t,t2 -> tuple(experiment, 'None')}.set{not_deconvoluted2}
        // file_cellmetadata = MULTIPLET.out.file__cellmetadata
        
        split_channel = vireo_out_sample_donor_ids.combine(ch_experiment_filth5, by: 0)
        split_channel2 = not_deconvoluted2.combine(ch_experiment_filth5, by: 0)
        // combining these 2 channels in one
        split_channel3 = split_channel.mix(split_channel2)
        // adding the scrublet paths to the channel.
        split_channel4 = split_channel3.combine(scrublet_paths, by: 0)

        split_channel5 = split_channel4.map{
            val_sample, val_donor_ids_tsv, val_filtered_matrix_h5, path_scrublet ->
            [  val_sample,file(val_donor_ids_tsv),file(val_filtered_matrix_h5),path_scrublet,params.outdir]
        }

        SPLIT_DONOR_H5AD(split_channel5)

        // collect file paths to h5ad files in tsv tables:
        header_seed="experiment_id\tdonor\th5ad_filepath"
        collect_file2(SPLIT_DONOR_H5AD.out.donors_h5ad_tsv.collect(),"donors_h5ad.tsv",params.outdir+'/deconvolution/filepaths',0,header_seed)

        header_seed="experiment_id\th5ad_filepath"
        collect_file3(SPLIT_DONOR_H5AD.out.exp__donors_h5ad_tsv.collect(),"exp__donors_h5ad.tsv",params.outdir+'/deconvolution/filepaths',0,header_seed)

        header_seed="experiment_id\tdonor\th5ad_filepath"
        collect_file4(SPLIT_DONOR_H5AD.out.donors_h5ad_assigned_tsv.collect(),"donors_h5ad_assigned.tsv",params.outdir+'/deconvolution/filepaths',0,header_seed)

        header_seed="experiment_id\th5ad_filepath"
        out_h5ad = collect_file5(SPLIT_DONOR_H5AD.out.exp__donors_h5ad_assigned_tsv.collect(),"exp__donors_h5ad_assigned.tsv",params.outdir+'/deconvolution/filepaths',0,header_seed)

        header_seed="experiment_id\th5ad_filepath"
        collect_file6(SPLIT_DONOR_H5AD.out.h5ad_tsv.collect(),"cellranger_as_h5ad.tsv",params.outdir+'/deconvolution/filepaths',0,header_seed)

        header_seed="experiment_id\tdonor\tn_cells"
        ch_vireo_donor_n_cells_tsv = collect_file7(REPLACE_GT_DONOR_ID2.out.sample_summary_tsv.collect(),"vireo_donor_n_cells.tsv",params.outdir+'/deconvolution/filepaths',0,header_seed)

        // paste experiment_id and donor ID columns with __ separator
        header_seed="experiment_id\tn_cells"
        vireo_out_sample__exp_summary_tsv = collect_file8(SPLIT_DONOR_H5AD.out.donor_n_cells.collect(),"vireo_exp__donor_n_cells.tsv",params.outdir+'/deconvolution/filepaths',0,header_seed)

        if (params.genotype_input.run_with_genotype_input & params.genotype_input.posterior_assignment) {
            if (params.extra_sample_metadata!='' & params.add_donor_metadata){
                // Here we have sample level metadata but we have chosen to keep the donor ids.
                // in this scenario we enhance the donor level vireo metadata file to add the donor metadata to the h5ads and eventually to the donor and teanche report
                ENHANCE_VIREO_METADATA_WITH_DONOR(params.extra_sample_metadata,vireo_out_sample__exp_summary_tsv,REPLACE_GT_DONOR_ID.out.assignments.collect())
                vireo_out_sample__exp_summary_tsv = ENHANCE_VIREO_METADATA_WITH_DONOR.out.replaced_vireo_exp__donor_n_cells_out
            }
        }


        PLOT_DONOR_CELLS(ch_vireo_donor_n_cells_tsv)



    emit:
        out_h5ad
        vireo_out_sample__exp_summary_tsv
        vireo_out_sample_donor_vcf
        ch_poolid_csv_donor_assignments = ch_poolid_donor_assignment
        sample_possorted_bam_vireo_donor_ids = ch_experiment_bam_vireo_donor_ids
        assignments_all_pools
}
