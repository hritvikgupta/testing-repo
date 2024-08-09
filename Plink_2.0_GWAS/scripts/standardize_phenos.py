import argparse as ap
import pandas as pd

from pathlib import Path
# from sklearn.preprocessing import StandardScaler


def make_arg_parser():
    parser = ap.ArgumentParser(description=".")

    parser.add_argument('-p', '--phenoCovarTable', required=True, help='.tsv Phenotype and covariate file. Should have IID column.')
    #parser.add_argument('-s', '--samples', required=True, help='.csv of cohort assignments')
    parser.add_argument('-s', '--samples', required=True, help='.tsv with sample IDs. Should have IID and FID columns.')
    parser.add_argument('-o','--outfile',default=None, help='Name of standardized output file')

    return parser


def get_basename(filepath, parent=False, suffixes=None):
    """
    Takes a path (string or PosixPath object) and returns the filename without path or suffix. 

    Args:
        filepath (filepath): filepath to operate on
        parent (bool, optional): Whether to keep parent path in name. Defaults to False.
        suffixes (list, optional): list of suffixes to strip (.txt, .csv), otherwise strips anything after last period. Defaults to None.

    Returns:
        str: stripped filename
    """
    filename = Path(filepath)
    # if suffixes supplied, only strip those
    if suffixes:
        while filename.suffix in set(suffixes):
            filename = filename.with_suffix("")
    else:
        filename = filename.with_suffix("")
    # if parent set to true, keep it
    if parent:
        filename = str(filename)
    else:
        filename = filename.name
    return filename

# parse arguments
args = make_arg_parser().parse_args()

samplefile = args.samples
pheno_covar_file = args.phenoCovarTable
outfile = args.outfile

# read in the datafiles
df = pd.read_table(pheno_covar_file, index_col=['FID', 'IID'])
samples = open(samplefile).read().splitlines()

# subsample by ID
df = df[df.index.get_level_values('IID').isin(samples)]
print(f"\nNumber of selected samples: {df.shape[0]}")

# let's categorize the column types
binary_columns = [col for col in df.columns if len(df[col].unique()) <= 3]
numerical_columns = df.select_dtypes(include=['float64', 'int64']).columns.to_list()
cat_columns = [col for col in df.columns if col not in numerical_columns]
cat_columns = cat_columns + binary_columns
quant_columns = [col for col in numerical_columns if col not in cat_columns]
print(f'Total Number of Columns: {df.shape[1]}\nNumber of category columns: {len(cat_columns)}\nNumber of quant columns: {len(quant_columns)}')

# Scale the data
# scaler = StandardScaler()

# I don't thinK we should scale the entire dataset since these are targets, not features...
# df[quant_columns] = scaler.fit_transform(df[quant_columns])

# instead Let's scale each column individually (Apply not working...)
# df[quant_columns] = df[quant_columns].apply(lambda x: StandardScaler().fit_transform(x))
# for col in quant_columns:
#     df[col] = scaler.fit_transform(df[[col]])

# let's scale using raw pandas
df[quant_columns] = (df[quant_columns] - df[quant_columns].mean()) / df[quant_columns].std()

df[binary_columns] = df[binary_columns].replace({0: 1, 1: 2}).fillna(-9)

# save
if outfile:
    df.reset_index().to_csv(outfile, sep='\t', index=False)
else:
    # outfile = base = get_basename(pheno_covar_file,parent=False) + '_standardized.tsv'
    outfile = 'plink2_pheno_covars_standardized.tsv'
    df.reset_index().to_csv(outfile, sep='\t', index=False, na_rep='NA')