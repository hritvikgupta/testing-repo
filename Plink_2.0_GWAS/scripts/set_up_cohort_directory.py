import pandas as pd
import argparse as ap

def make_arg_parser():
    parser = ap.ArgumentParser(description=".")
    
    parser.add_argument('-d', '--data', required=True, help='.csv Phenotype and covariate file')
    parser.add_argument('-c', '--cohort', required=True, help='Cohort to set up')
    parser.add_argument('-s', '--samples', required=True, help='.csv of cohort assignments')
    parser.add_argument('-i', '--id', required=True, help='Column with sample IDs')
    parser.add_argument('--plinkFam', required=True)

    return parser

args = make_arg_parser().parse_args()

id_col = args.id
cohort = args.cohort
plink_fam = args.plinkFam

data = pd.read_csv(args.data, index_col=id_col, dtype={id_col: str})
samples = pd.read_csv(args.samples, index_col=id_col, dtype={id_col: str})

print(data)
print(samples)

plink_fam = pd.read_table(plink_fam, header=None, comment='#', index_col=1, sep='\\s+', dtype=str)

keep_samples = data.index.intersection(samples.index).intersection(plink_fam.index)

data, samples = data.loc[keep_samples], samples.loc[keep_samples]
cohort_samples = samples.index[samples[cohort] == 1]

open('sample_list.txt', 'w+').write('\n'.join(cohort_samples))

data = data.loc[cohort_samples]
data.index.name = 'FID'
data.insert(0, 'IID', data.index)
data.to_csv('plink2_pheno_covars.txt', sep='\t')
