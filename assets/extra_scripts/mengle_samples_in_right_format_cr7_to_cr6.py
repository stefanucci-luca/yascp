import glob
import pandas as pd
import sys
import os
# This script mengles the data in the right format for use in YASCP pipeline.
# this converts the cellranger cellranger7.10 format to cellranger6.1.1

Cellranger_path = '/lustre/scratch123/hgi/mdt1/projects/immunodeficiency/trego/cellranger_output/cellranger_data'
# Cellranger_path = '/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/ania/test_data'
Outdir = '/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/ania/mengled_samples'
all_files = glob.glob(f'{Cellranger_path}/cellranger710*')
Sample_names = pd.read_csv('/lustre/scratch123/hgi/mdt1/projects/immunodeficiency/trego/cellranger_output/genotype_deconvolution/sample_metadata.txt',sep=' ')
# Sample_names = Sample_names.set_index('file')
Dataset3 = {}
c1 = 0
for file1 in all_files:
    print(file1)
    c1 +=1
    
    Sample_Name = file1.split('/')[-1].split('_')[-1]
    Sample_name2 = ' or '.join(Sample_names[Sample_names['file'].str.contains(Sample_Name)]['sample'])
    Sampe_multiplex = Sample_names[Sample_names['file'].str.contains(Sample_Name)]
    os.mkdir(f'{Outdir}/{Sample_name2}')
    Metrics_file = f"{file1}/per_sample_outs/cellranger710_multi_{Sample_Name}/metrics_summary.csv"
    Dataset = pd.read_csv(Metrics_file)
    Dataset1 = Dataset[Dataset['Library Type']=='Gene Expression']
    Dataset1 = Dataset1.set_index('Metric Name')
    Dataset1 = Dataset1['Metric Value']
    Dataset2 = pd.DataFrame(Dataset1).T
    Dataset2.to_csv(f'{Outdir}/{Sample_name2}/metrics_summary.csv',index=False)
    'per_sample_outs/cellranger710_multi_05d72fb896018b97366286cde31ce7d6/count/sample_alignments.bam'
    os.symlink(f"{file1}/per_sample_outs/cellranger710_multi_{Sample_Name}/count/analysis", f"{Outdir}/{Sample_name2}/analysis")
    os.symlink(f"{file1}/per_sample_outs/cellranger710_multi_{Sample_Name}/count/sample_filtered_feature_bc_matrix", f"{Outdir}/{Sample_name2}/filtered_feature_bc_matrix")
    os.symlink(f"{file1}/multi/count/raw_feature_bc_matrix", f"{Outdir}/{Sample_name2}/raw_feature_bc_matrix")
    os.symlink(f"{file1}/per_sample_outs/cellranger710_multi_{Sample_Name}/count/sample_filtered_feature_bc_matrix.h5", f"{Outdir}/{Sample_name2}/filtered_feature_bc_matrix.h5")
    os.symlink(f"{file1}/per_sample_outs/cellranger710_multi_{Sample_Name}/count/sample_molecule_info.h5", f"{Outdir}/{Sample_name2}/molecule_info.h5")
    os.symlink(f"{file1}/per_sample_outs/cellranger710_multi_{Sample_Name}/count/sample_alignments.bam", f"{Outdir}/{Sample_name2}/possorted_genome_bam.bam")
    os.symlink(f"{file1}/per_sample_outs/cellranger710_multi_{Sample_Name}/count/sample_alignments.bam.bai", f"{Outdir}/{Sample_name2}/possorted_genome_bam.bam.bai")
    os.symlink(f"{file1}/multi/count/raw_feature_bc_matrix.h5", f"{Outdir}/{Sample_name2}/raw_feature_bc_matrix.h5")
    os.symlink(f"{file1}/per_sample_outs/cellranger710_multi_{Sample_Name}/web_summary.html", f"{Outdir}/{Sample_name2}/web_summary.html")
    Dataset3[c1]={'experiment_id':Sample_name2,'n_pooled':1+Sampe_multiplex['spikein'].values[0]+Sampe_multiplex['chimeric'].values[0],'donor_vcf_ids':"''",'data_path_10x_format':f"{Outdir}/{Sample_name2}"}
Yascp_input = pd.DataFrame(Dataset3).T  
Yascp_input.to_csv(f'{Outdir}/input.tsv',sep='\t',index=False)
print('Done')


