#!/usr/bin/env python

__date__ = '2021-07-28'
__version__ = '0.0.1'

import os
from shutil import copyfile,copytree,copy
from os import listdir
import glob
import argparse
import pandas as pd



def main_data_colection(pipeline='',name='',directory='',input_table=None,cb_res=None,web_transfer=False,project_name='all'):
    
    name_dir=f"Summary_plots/{name}"
    try:
        os.makedirs(f'{name_dir}')
    except:
        print("dir exists")


    #Cellbender
    print('prepearing Cellbender folder')
    try:
        os.mkdir(f'{name_dir}/Cellbender')
    except:
        print('dire exists')
    
    if (cellbender)=='cellranger':
        # here we do not use cellbender and go with default cellranger
        df_cellbender = None
    elif (cellbender)=='cellbender':
        # here we have run the cellbender as par of pipeline. 
        # file_path = glob.glob(f'{results_dir}/nf-preprocessing/cellbender/qc_cluster_input_files/*{cb_res}*')[0]
        file_path = glob.glob(f'{results_dir}/nf-preprocessing/cellbender/*/cellbender-epochs_*/*{cb_res}*10x_mtx*')
        file_path2 = glob.glob(f'{results_dir}/nf-preprocessing/cellbender/*/*{cb_res}*10x_mtx*')
        joined_file_paths = file_path+file_path2
        df_cellbender = pd.DataFrame(joined_file_paths,columns=['data_path_10x_format'])
        df_cellbender['experiment_id']=df_cellbender['data_path_10x_format'].str.split('/').str[3]
        df_cellbender= df_cellbender.set_index('experiment_id')
    
    print(df_cellbender)
    # folder1 = f'{directory}/nf-preprocessing/cellbender'
    # folder2 = f'{directory}/cellbender_vs_cellranger'
    
    resolution2=cb_res.replace('pt','.')
    if (df_cellbender is not None):
        # print('t')
        for folder in df_cellbender.index:
            print(folder)
            print("yes!!")
            # dir1 = f"{df_cellbender.loc[folder,'data_path_10x_format']}/.."
            dir1 = '/'.join(df_cellbender.loc[folder,'data_path_10x_format'].split('/')[:-1])
            dir = '/'.join(df_cellbender.loc[folder,'data_path_10x_format'].split('/')[:-2])
            print(dir1)
            if os.path.isdir(dir1):
                print("yes22!!")
                try:
                    copyfile(f'{dir1_b}/plots/cellbender_results-cellbender_FPR_{cb_res}_filtered-ambient_signature-scatter_genenames.png', f'{name_dir}/Cellbender/{folder}_ambient_signature-scatter_genenames.png')
                except:
                    print('missing1')
                try:
                    copyfile(f'{dir1}/plots/cellbender_results-cellbender_FPR_{cb_res}_filtered-abs_count_difference-boxplot.png', f'{name_dir}/Cellbender/{folder}_ount_difference-boxplot.png')
                except:
                    print('missing2')
                try:
                    copyfile(f'{dir1}/plots/cellbender.pdf', f'{name_dir}/Cellbender/cellbender_{folder}.pdf')
                except:
                    print('missing3')
                try:
                    copyfile(f'{dir}/compare_cellranger/fpr_{resolution2}/boxplots_cellranger_vs_cellbender.png', f'{name_dir}/Cellbender/{folder}_boxplots_cellranger_vs_cellbender.png')
                except:
                    print('missing4')                
                try:
                    copyfile(f'{dir}/compare_cellranger/fpr_{resolution2}/barcode_vs_total_counts.png', f'{name_dir}/Cellbender/{folder}_barcode_vs_total_counts.png')
                except:
                    print('missing5')                
                try:
                    copyfile(f'{dir}/compare_cellranger/fpr_{resolution2}/boxplot_topgenes_cellranger_vs_cellbender.png', f'{name_dir}/Cellbender/{folder}_boxplot_topgenes_cellranger_vs_cellbender.png')
                except:
                    print('missing6')
    # Fetch Gather
    df_raw = pd.read_table(input_table, index_col = 'experiment_id')
    # 'Note that the names for the future projects may be different - have to be handled on the Nextflow modules'
    print('prepearing fetch folder')
    try:
        os.mkdir(f'{name_dir}/Fetch Pipeline')
    except:
        print('exists')
    metadata_table = pd.DataFrame()
    d1 = glob.glob(f"{results_dir}/*/web_summary.html")
    d2 = glob.glob(f"{results_dir}/*/*/web_summary.html")
    d3 = glob.glob(f"{results_dir}/*/*/*/web_summary.html")
    d4 =glob.glob(f"{results_dir}/*/*/*/*/web_summary.html")
    ts3 = [*d1, *d2,*d3,*d4]
    ts3 = pd.DataFrame(ts3,columns=['col'])
    for folder in df_raw.index:
        dir3 = f"{df_raw.loc[folder, 'data_path_10x_format']}"
        try:
            copyfile(f'{dir3}/web_summary.html', f'{name_dir}/Fetch Pipeline/html_{folder}.html')
            metadata = pd.read_csv(f'{dir3}/metrics_summary.csv',sep=',',index_col=False)
            metadata['Sample_id']=folder
            metadata.set_index('Sample_id',drop=True,inplace=True)
            metadata_table=pd.concat([metadata_table,metadata])
        except:
            try:
                try:
                    copyfile(f'{results_dir}/handover/{name_dir}/Fetch Pipeline/html_{folder}.html', f'{name_dir}/Fetch Pipeline/html_{folder}.html')
                except:
                    try:
                        copyfile(ts3[ts3['col'].str.contains(f'{folder}')]['col'].values[0], f'{name_dir}/Fetch Pipeline/html_{folder}.html')
                    except:    
                        copyfile(f'/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/{name}/results_rsync2/results/handover/{name_dir}/Fetch Pipeline/html_{folder}.html', f'{name_dir}/Fetch Pipeline/html_{folder}.html')
            except:
                _='pass'
            try:
                metadata_table = pd.read_csv(f'{results_dir}/handover/{name_dir}/Fetch Pipeline/CSV/Submission_Data_Pilot_UKB.file_metadata.tsv',sep='\t')
                metadata_table.set_index('Sample_id',drop=True,inplace=True)
            except:
                try:
                    metadata_tablemetadata_table = pd.read_csv(f'{results_dir}/handover/{name_dir}/Fetch Pipeline/Submission_Data_Pilot_UKB.file_metadata.tsv',sep=',')
                    metadata_table.set_index('Sample_id',drop=True,inplace=True)
                except:
                    metadata_table = pd.DataFrame()
    try:
        os.mkdir(f'{name_dir}/Fetch Pipeline/CSV')
    except:
        print('exists')
    if (len(metadata_table)>0):
        metadata_table.to_csv(f'{name_dir}/Fetch Pipeline/CSV/Submission_Data_Pilot_UKB.file_metadata.tsv',sep='\t')
    else:
        try:
            copyfile(f'/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/{name}/results_rsync2/results/handover/Summary_plots/{name}/Fetch Pipeline/CSV/Submission_Data_Pilot_UKB.file_metadata.tsv', f'{name_dir}/Fetch Pipeline/CSV/Submission_Data_Pilot_UKB.file_metadata.tsv')
        except:
            _='pass' 

    try:
        os.mkdir(f'{name_dir}/Fetch Pipeline/Input')
        copyfile(input_table, f'{name_dir}/Fetch Pipeline/Input/input_table.tsv')
    except:
        print('exists')

    try:
        os.mkdir(f'{name_dir}/GT Match___1000')
    except:
        print('exists')

    #NOW COPY THE GT OUTPUTS 
    #folder1 = f'{directory}/deconvolution/vireo_gt_fix'
    #if os.path.isdir(folder1):
    #    copyfile(f'{folder1}/assignments_all_pools.tsv', f'{name_dir}/GT Match___1000/assignments_all_pools.tsv')
        

    #if (pipeline=='Deconvolution'):
    folder1 = f'{directory}/deconvolution/split_donor_h5ad'
    if os.path.isdir(folder1):
        print('prepearing Deconvolution folder')
        try:
            os.mkdir(f'{name_dir}/Deconvolution')
        except:
            print('dire exists')
        Folders = listdir(folder1)
        for folder in Folders:
            copyfile(f'{folder1}/{folder}/Vireo_plots.pdf', f'{name_dir}/Deconvolution/Vireo_plots_{folder}.pdf')
        
    folder1 = f'{directory}/celltype/celltypist'
    if os.path.isdir(folder1):
        print('prepearing celltype folder')
        try:
            os.mkdir(f'{name_dir}/Cell-type assignment')
        except:
            print('dire exists')
        try:
            os.mkdir(f'{name_dir}/Cell-type assignment/celltypist')
        except:
            print('dire exists')
        Folders = listdir(folder1)
        for model_type in Folders:
            print(model_type)
            Folders2 = listdir(f"{folder1}/{model_type}")
            for donor in Folders2:
                copyfile(f'{folder1}/{model_type}/{donor}/{donor}_predicted_labels.pdf', f'{name_dir}/Cell-type assignment/celltypist/{model_type}_{donor}_predicted_labels.pdf')
                copyfile(f'{folder1}/{model_type}/{donor}/{donor}_majority_voting.pdf', f'{name_dir}/Cell-type assignment/celltypist/{model_type}_{donor}_majority_voting.pdf')

    folder1 = f'{directory}/celltype/azimuth'


    # folder1 = f'{directory}/plots'
    if os.path.isdir(folder1):
        try:
            os.mkdir(f'{name_dir}/Cell-type assignment')
        except:
            print('dire exists')
        # try:
        #     os.mkdir(f'{name_dir}/Cell-type assignment/azimuth')
        # except:
        #     print('dire exists')
        # copyfile(fil1, f'{name_dir}/QC metrics/plot_ecdf-x_log10.var=total_counts.color=experiment_id-adata.png')
        files = glob.glob(f'{folder1}/*[!.gz]')
        for file1 in files:
            # print(file1)
            copy(file1, f'{name_dir}/Cell-type assignment')


    folder1 = f'{directory}/plots/per_celltype_outliers'
    if os.path.isdir(folder1):
        print('yes')
        try:
            os.mkdir(f'{name_dir}/QC metrics')
        except:
            print('dire exists')
        os.system(f'ln -s ./{folder1} {name_dir}/QC metrics')


    folder1 = f'{directory}/plots'
    if os.path.isdir(folder1):
        try:
            os.mkdir(f'{name_dir}/QC metrics')
        except:
            print('dire exists')
            
        density_files = glob.glob(f'{folder1}/*cell_desity*')
        density_files2 = glob.glob(f'{folder1}/merged_h5ad/plots/*cell_desity*')
        density_files.extend(density_files2)
        for dens_file in density_files:
            try:
                os.mkdir(f'{name_dir}/QC metrics/density')
            except:
                print('dire exists')
            nam1 = dens_file.split('/')[-1]
            copyfile(f'{dens_file}', f'{name_dir}/QC metrics/density/{nam1}')

        density_files = glob.glob(f'{folder1}/*adata-outlier_cells*')
        density_files2 = glob.glob(f'{folder1}/merged_h5ad/plots/*adata-outlier_cells*')
        density_files.extend(density_files2)
        for dens_file in density_files:
            try:
                os.mkdir(f'{name_dir}/QC metrics/outlier_cells')
            except:
                print('dire exists')
            nam1 = dens_file.split('/')[-1]
            copyfile(f'{dens_file}', f'{name_dir}/QC metrics/outlier_cells/{nam1}')
        
        density_files = glob.glob(f'{folder1}/*per_celltype_outliers*')
        for dens_file in density_files:
            print(dens_file)
            name_dir_final = dens_file.split('/')[-1]
            copytree(dens_file, f'{name_dir}/QC metrics/{name_dir_final}')         
        copyfile(f'{folder1}/adata-cell_filtered_per_experiment-n_cells_before_after.png', f'{name_dir}/QC metrics/adata-cell_filtered_per_experiment-n_cells_before_after.png')
        copyfile(f'{folder1}/scatterplot-sex_sample_swap_check.png', f'{name_dir}/QC metrics/scatterplot-sex_sample_swap_check.png')
        fil1 = glob.glob(f'{folder1}/plot_ecdf-x_log10*total_counts*')[0]
        copyfile(fil1, f'{name_dir}/QC metrics/plot_ecdf-x_log10.var=total_counts.color=experiment_id-adata.png')


    folder1 = f'{directory}/clustering'
    if os.path.isdir(folder1):
        try:
            os.mkdir(f'{name_dir}/Clustering')
            os.mkdir(f'{name_dir}/Clustering/Harmony')
        except:
            print('dire exists')
            try:
                os.mkdir(f'{name_dir}/Clustering/Harmony')
            except:
                print('dire exists')
        Harmony_UMAPS = glob.glob(f'{folder1}/*/*harmony*/*/plots/umap*')
        for umap1 in Harmony_UMAPS:
            
            name = umap1.split('/')[-1]
            resolution = umap1.split('/')[-3].split("resolution=")[1]
            
            copyfile(umap1, f'{name_dir}/Clustering/Harmony/{resolution}res_Harmony_{name}')


    folder1 = f'{directory}/clustering'
    if os.path.isdir(folder1):
        try:
            os.mkdir(f'{name_dir}/Clustering')
            os.mkdir(f'{name_dir}/Clustering/BBKNN')
        except:
            print('dire exists')
            try:
                os.mkdir(f'{name_dir}/Clustering/BBKNN')
            except:
                print('dire exists')
        Harmony_UMAPS = glob.glob(f'{folder1}/*/*bbknn*/*/plots/umap*')
        for umap1 in Harmony_UMAPS:
            
            name = umap1.split('/')[-1]
            resolution = umap1.split('/')[-3].split("resolution=")[1]
            
            copyfile(umap1, f'{name_dir}/Clustering/BBKNN/{resolution}res_BBKNN_{name}')


    folder1 = f'{directory}/UMAPs'
    if os.path.isdir(folder1):
        try:
            os.mkdir(f'{name_dir}/Clustering')
            os.mkdir(f'{name_dir}/Clustering/Coloured')
        except:
            try:
                os.mkdir(f'{name_dir}/Clustering/Coloured')
            except:
                print('dire exists')
            print('dire exists')
        Coloured_UMAPS = glob.glob(f'{folder1}/*')
        for umap1 in Coloured_UMAPS:
            name = umap1.split('/')[-1]
            copyfile(umap1, f'{name_dir}/Clustering/Coloured/{name}')                    
    
    
    folder1 = f'{directory}/handover/Donor_Quantification_summary'
    if os.path.isdir(folder1):
        try:
            copytree(folder1, f'{name_dir}/Summary')
        except:
            print('dire exists')
    folder1 = f'{directory}/celltype'
    if os.path.isdir(folder1):
        try:
            copyfile(f"{folder1}/donor_celltype_report.tsv", f'{name_dir}/Summary/donor_celltype_report.tsv')  
            copyfile(f"{folder1}/tranche_celltype_report.tsv", f'{name_dir}/Summary/tranche_celltype_report.tsv') 
        except:
            _ ="Cant copy"
        
    if web_transfer !='false':
        print('Here we triger the CI placed script, which happens only for the local Cardinal project')
        os.system(f'exit && bash ../../../scripts/rsync_to_web.sh {project_name}')

    # if secret params exist then transfer data to web. 
    # rsync -vr Summary_plots ubuntu@1xxxxx:/volume/scRNA_test_app/scrna_static_and_media_files/media

if __name__ == "__main__":
    # This code will gather the most important graphs and place them in a summary folder.
    parser = argparse.ArgumentParser(
        description="""
            Collect important graphs for web visualisations
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-rd', '--results_dir',
        action='store',
        dest='results_dir',
        required=True,
        help='The nextflow results directory.'
    )

    parser.add_argument(
        '-pn', '--project_name',
        action='store',
        dest='project_name',default='all',
        required=False,
        help='Project name to use for web transfer.'
    )

    parser.add_argument(
        '-cb_res', '--cb_res',
        action='store',
        dest='cb_res',
        required=True,
        help='The cellbender resolution used in downstream analysis'
    )

    parser.add_argument(
        '-wt', '--web_transfer',
        action='store',
        dest='web_transfer',
        required=False,default=False,
        help='If run internally we will transfer the results to a website if run through Gitlab CI'
    )

    parser.add_argument(
        '-it', '--input_table',
        action='store',
        dest='input_table',
        required=True,
        help='Input table with paths to 10xRuns'
    )

    parser.add_argument(
        '-cb', '--cellbender',
        action='store',
        dest='cellbender',
        required=True,
        help='The path to cellbender runs, since this may be external.'
    )
    options = parser.parse_args()
    
    t=os.getcwd()
    results_dir = options.results_dir
    name = t.split('/')[-4]
    
    input_table=options.input_table
    cb_res=options.cb_res
    cellbender = options.cellbender
    web_transfer = options.web_transfer
    project_name = options.project_name
    main_data_colection(name=f"{name}",directory=results_dir,input_table=input_table,cb_res=cb_res,web_transfer=web_transfer,project_name=project_name)
