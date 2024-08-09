from manhattan_plot import ManhattanPlot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse as ap
import os


def make_arg_parser():
    parser = ap.ArgumentParser(description=".")

    # Add a non-optional argument for phenotype
    parser.add_argument('-p', '--phenotype', required=True, help='phenotype')
    # Add a non-optional argument for cohort
    parser.add_argument('-c', '--cohort', required=True,
                        help='which cohort')
    # Add an argument for output_directory
    parser.add_argument('-o', '--outDir', default='./',
                        help='Path to output directory. Default: current working directory')
    # Add an argument for summary statistics file
    parser.add_argument('-s', '--sumstats', required=True,
                        help='Path to summary statistics file')
    # add argument for addint annotationgs
    parser.add_argument('-a', '--annot', required=False, default=None)
    parser.add_argument('-col', '--colnames', required=True, help='File with column name mappings')


    return parser

def make_column_name_map_dict(colnames_file):
    colnames_rows = open(colnames_file).read().splitlines()
    print(colnames_rows)
    col_map = dict(zip([r.split('=')[0] for r in colnames_rows],
                    [r.split('=')[1] for r in colnames_rows]))
    print(col_map)
    return(col_map)

# parse arguments
args = make_arg_parser().parse_args()
output_dir = args.outDir
pheno, cohort = args.phenotype, args.cohort
# pheno_table = args.phenoTable # only used for trait cardinality for plot
sumstats_file = args.sumstats
annot_file = args.annot


output_manhattan = f'{output_dir}/{cohort}.{pheno}.manhattan.png'
output_qq = f'{output_dir}/{cohort}.{pheno}.qq.png'

# Example output file to be used lives in:
# /project/pmbb_codeworks/projects/geno_pheno_workbench_dev/GWAMA_META/Meta/Sumstats/*.gz

# Instantiate manhattan plot object
plot_title = f'Plink2 GWAS Manhattan for {cohort}: {pheno.replace("_", " ")}'
mp = ManhattanPlot(sumstats_file, title=plot_title)
mp.load_data()
# clean data, use parameter map function
plink2_col_map = make_column_name_map_dict(args.colnames)
neat_col_map = {plink2_col_map['#CHROM']: '#CHROM',
                plink2_col_map['POS']: 'POS',
                plink2_col_map['ID']: 'ID',
                plink2_col_map['P']: 'P'}

mp.clean_data(col_map=neat_col_map)

# add conditional adventure for annotations
if annot_file is not None:
    annot_df = pd.read_csv(annot_file)
    annot_df['ID'] = annot_df['Gene']
    mp.add_annotations(annot_df, extra_cols=['RSID'])

mp.get_thinned_data()
# mp.thinned = mp.thinned.dropna(subset='P')
if ~np.any(mp.thinned['P'] < 5E-8):
    p_thresh = np.nanquantile(mp.thinned['P'], 10 / len(mp.thinned))
else:
		p_thresh = 5E-8

mp.update_plotting_parameters(vertical=True,sig=p_thresh,sug=p_thresh,annot_thresh=p_thresh,merge_genes=True)


# mp.update_plotting_parameters(vertical=True, merge_genes=True)
# mp.full_plot(save=output_manhattan,rep_genes=known_genes,rep_boost=True)
mp.full_plot(save=output_manhattan,rep_boost=True)

# close fig and save
plt.clf()
print(f"Saved Manhattan plot to: {output_manhattan}")

# mp.qq_plot
mp.qq_plot(save=output_qq)
print(f"Saved qq plot to: {output_qq}")
