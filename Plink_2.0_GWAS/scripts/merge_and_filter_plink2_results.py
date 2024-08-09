import pandas as pd
import argparse as ap

def make_arg_parser():
    parser = ap.ArgumentParser(description=".")
    
    # Add non-optional list arguments for phenotypes
    parser.add_argument('-s', '--sumstats', nargs='+', help='List of results files')
    parser.add_argument('-c', '--colnames', required=True, help='File with column name mappings')
    parser.add_argument('-p', '--pheno', help='Phenotype')
    parser.add_argument('--cohort', help='Cohort')
    parser.add_argument('--pvalue', help='P-value for filtering', type=float, default=1E-5)
    
    return parser

args = make_arg_parser().parse_args()
colnames_file = args.colnames
input_files = args.sumstats
p_thresh = args.pvalue

pheno = args.pheno
cohort = args.cohort

merge_output = f'{pheno}.plink2.gz'
filter_output = f'{pheno}.filtered.plink2.csv'

colnames_rows = open(colnames_file).read().splitlines()
print(colnames_rows)
col_map = dict(zip([r.split('=')[0] for r in colnames_rows],
                   [r.split('=')[1] for r in colnames_rows]))

print(col_map)

dfs = []

for f in input_files:
    print(f)
    f_pheno = '.'.join(f.split('/')[-1].split('.')[:-3])
    temp = pd.read_table(f)
    temp['PHENO'] = f_pheno
    dfs.append(temp)

all = pd.concat(dfs)
print(all)
all['A2'] = all.apply(lambda x: x['REF'] if x['A1'] == x['ALT'] else x['ALT'], axis=1)
print(all)
print(all.columns)

all = all.rename(columns=col_map)
all.to_csv(merge_output, sep='\t', index=False, na_rep='NA')


all_filtered = all[all[col_map['P']] <= p_thresh]
all_filtered['COHORT'] = cohort
all_filtered['PHENO'] = pheno
all_filtered.to_csv(filter_output, index=False, na_rep='NA')