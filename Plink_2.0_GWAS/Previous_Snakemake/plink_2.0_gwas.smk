include: '/project/ritchie/lindsay_snakemake_workflows/liftover_wrapper/Snakefile'
include: '/project/ritchie/lindsay_snakemake_workflows/biofilter_wrapper/Snakefile'

configfile: 'config_plink_2.0_gwas.yaml'

AUTOSOMES=list(range(1, 23))

wildcard_constraints:
    cohort_dir='(' + '|'.join(config['plink2_cohorts']) + ')',
    pheno='(' + '|'.join(config['plink_bin_phenos']) + '|' + '|'.join(config['plink_quant_phenos']) +')'

rule standardize_phenos:
    output:
        '{cohort_dir}/plink2_{data}_standardized.txt'
    input:
        data='{cohort_dir}/plink2_{data}.txt',
        samples='{cohort_dir}/sample_list.txt'
    run:
        import pandas as pd

        df = pd.read_csv(input.data, sep='\t', index_col=['FID', 'IID'])
        samples = open(input.samples).read().splitlines()

        df = df[df.index.get_level_values('IID').isin(samples)]

        non_bin_cols = [c for c in df.columns if len(df[c].unique()) > 2]
        print(df)
        print(non_bin_cols)

        df[non_bin_cols] = df[non_bin_cols].apply(lambda x: (x - x.mean()) / x.std())

        df.to_csv(str(output), sep='\t', na_rep='NA')

rule call_plink2_logistic:
    output:
        '{cohort_dir}/GWAS_Results/bin_phenos.{chr}.{pheno}.glm.logistic.hybrid'
    input:
        pheno='{cohort_dir}/plink2_pheno_standardized.txt',
        covar='{cohort_dir}/plink2_covars_standardized.txt',
        plink_set=expand(config['plink_prefix'] + '{{chr}}{ext}', ext=config['plink_extensions']),
        samples='{cohort_dir}/sample_list.txt',
        longrange='/project/ritchie/datasets/Longrange_Regions/longrange_LD_hg' + config['build'] + '.bed'
    params:
        covar_list=lambda wildcards: ','.join([c for c in config['plink_covars'] if (wildcards.cohort_dir not in config['sex_strat_cohorts']) or (c not in config['sex_dependent_covars'])]),
        output_prefix='{cohort_dir}/GWAS_Results/bin_phenos.{chr}',
        plink_prefix=config['plink_prefix'] + '{chr}',
        plink_flag=config['plink_flag'],
        maf=config['plink_maf'],
        geno=config['plink_geno'],
        mind=config['plink_mind'],
        hwe=config['plink_hwe']
    envmodules:
        'plink/2.0-20210505'
    resources: mem_mb=16000
    threads: 4
    shell:
        """
        plink --glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+firth \
          --exclude range {input.longrange} \
          --keep-fam {input.samples} \
          --maf {params.maf} \
          --geno {params.geno} \
          --mind {params.mind} \
          --hwe {params.hwe} \
          --ci 0.95 \
          {params.plink_flag} {params.plink_prefix} \
          --pheno {input.pheno} \
          --pheno-name {wildcards.pheno} \
          --covar {input.covar} \
          --covar-name {params.covar_list} \
          --out {params.output_prefix}
        """

rule call_plink2_linear:
    output:
        '{cohort_dir}/GWAS_Results/quant_phenos.{chr}.{pheno}.glm.linear'
    input:
        pheno='{cohort_dir}/plink2_pheno_standardized.txt',
        covar='{cohort_dir}/plink2_covars_standardized.txt',
        plink_set=expand(config['plink_prefix'] + '{{chr}}{ext}', ext=config['plink_extensions']),
        samples='{cohort_dir}/sample_list.txt',
        longrange='/project/ritchie/datasets/Longrange_Regions/longrange_LD_hg' + config['build'] + '.bed'
    params:
        covar_list=lambda wildcards: ','.join([c for c in config['plink_covars'] if (wildcards.cohort_dir not in config['sex_strat_cohorts']) or (c not in config['sex_dependent_covars'])]),
        output_prefix='{cohort_dir}/GWAS_Results/quant_phenos.{chr}',
        plink_prefix=config['plink_prefix'] + '{chr}',
        plink_flag=config['plink_flag'],
        maf=config['plink_maf'],
        geno=config['plink_geno'],
        mind=config['plink_mind'],
        hwe=config['plink_hwe']
    envmodules:
        'plink/2.0-20210505'
    resources: mem_mb=16000
    threads: 4
    shell:
        """
        plink --glm hide-covar cols=+a1freq \
          --exclude range {input.longrange} \
          --keep-fam {input.samples} \
          --maf {params.maf} \
          --geno {params.geno} \
          --mind {params.mind} \
          --hwe {params.hwe} \
          --ci 0.95 \
          {params.plink_flag} {params.plink_prefix} \
          --pheno {input.pheno} \
          --pheno-name {wildcards.pheno} \
          --covar {input.covar} \
          --covar-name {params.covar_list} \
          --out {params.output_prefix}
        """

def get_chr_merge_input(wildcards):
    quant_expand = expand('{{cohort_dir}}/GWAS_Results/quant_phenos.{chr}.{{pheno}}.glm.linear', chr=config['plink_chromosomes'])
    bin_expand = expand('{{cohort_dir}}/GWAS_Results/bin_phenos.{chr}.{{pheno}}.glm.logistic.hybrid', chr=config['plink_chromosomes'])
    return quant_expand if wildcards.pheno in config['plink_quant_phenos'] else bin_expand

rule merge_plink2_output:
    output:
        '{cohort_dir}/Sumstats/{pheno}.sumstats.gz'
    input:
        get_chr_merge_input
    resources: mem_mb=30000
    run:
        import pandas as pd

        dfs = []
        for f in input:
            print(f)
            dfs.append(pd.read_table(f, sep='\s+'))

        df = pd.concat(dfs)
        df = df.dropna(subset=['P'])
        print(df)
        df.to_csv(str(output), sep='\t', index=False)

rule get_plink2_suggestive_snps:
    output:
        '{cohort_dir}/Suggestive/{pheno}.txt'
    input:
        '{cohort_dir}/Sumstats/{pheno}.sumstats.gz'
    params:
        p_col = lambda wildcards: 15 if wildcards.pheno in config['plink_quant_phenos'] else 16
    shell:
        """
        zcat {input} | awk '{{if (NR==1 || ${params.p_col} <= 1E-5) {{print $1, $3, $2, ${params.p_col}}}}}' > {output}
        """

rule get_plink2_suggestive_snps_lifted:
    output:
        '{cohort_dir}/Suggestive_Lifted/{pheno}.txt'
    input:
        '{cohort_dir}/Sumstats/{pheno}.liftover.tsv.gz'
    params:
        p_col = lambda wildcards: 15 if wildcards.pheno in config['plink_quant_phenos'] else 16
    shell:
        """
        zcat {input} | awk '{{if (NR==1 || ${params.p_col} <= 1E-5) {{print $1, $2, $3, ${params.p_col}}}}}' > {output}
        """

rule get_plink2_suggestive_input:
    output:
        'Annotations/plink2_suggestive_biofilter_input_positions.txt'
    input:
        expand('{cohort_dir}/Suggestive/{pheno}.txt', cohort_dir=config['plink2_cohorts'], pheno=config['plink_bin_phenos']),
        expand('{cohort_dir}/Suggestive/{pheno}.txt', cohort_dir=config['plink2_cohorts'], pheno=config['plink_quant_phenos'])
    resources: mem_mb=12000
    shell:
        """
        cat {input} | awk '{{print $1, $2, $3}}' | sort -n | uniq > {output}
        """

rule get_plink2_suggestive_lifted_input:
    output:
        'Annotations/lifted_plink2_suggestive_biofilter_input_positions.txt'
    input:
        expand('{cohort_dir}/Suggestive_Lifted/{pheno}.txt', cohort_dir=config['plink2_cohorts'], pheno=config['plink_bin_phenos']),
        expand('{cohort_dir}/Suggestive_Lifted/{pheno}.txt', cohort_dir=config['plink2_cohorts'], pheno=config['plink_quant_phenos'])
    resources: mem_mb=12000
    shell:
        """
        cat {input} | awk '{{print $1, $3, $2}}' | sort -n | uniq > {output}
        """

rule make_plink2_sumstats_liftover_input:
    output:
        '{cohort_dir}/Sumstats/{pheno}.liftover_input.txt'
    input:
        sumstats='{cohort_dir}/Sumstats/{pheno}.sumstats.gz'
    shell:
        """
        zcat {input} | tail -n +2 | awk '{{print "chr"$1, $3, $3+1}}' > {output}
        """

rule join_plink2_liftover_with_sumstats:
    output:
        '{cohort_dir}/Sumstats/{pheno}.liftover.tsv.gz'
    input:
        lo_out='{cohort_dir}/Sumstats/{pheno}.liftover_output.txt',
        lo_failed='{cohort_dir}/Sumstats/{pheno}.liftover_failed.txt',
        table='{cohort_dir}/Sumstats/{pheno}.sumstats.gz'
    params:
        chrCol=1,
        posCol=3,
        has_header=True,
        keep_failed_rows=True
    resources:
        mem_mb=35000
    run:
        import pandas as pd
        import os

        ss = pd.read_table(input['table'], nrows=None, sep='\s+', header=None if not params.has_header else 'infer')
        chrColIndex = int(params['chrCol']) - 1
        posColIndex = int(params['posCol']) - 1
        original_column_order = ss.columns
        chrCol = ss.columns[chrColIndex]
        posCol = ss.columns[posColIndex]
        ss[chrCol] = ss[chrCol].astype(str)
        ss[posCol] = ss[posCol].astype(int)
        ss = ss.set_index([chrCol, posCol])
        print(ss)

        if os.path.getsize(str(input['lo_failed'])) != 0:
            failed = pd.read_table(input['lo_failed'], comment='#', header=None, names=['chrCHR', 'POS', 'END'])
            failed['CHR'] = failed['chrCHR'].str[3:]
        else:
            failed = pd.DataFrame(columns=['chrCHR', 'POS', 'END', 'CHR'])

        new = pd.read_table(input['lo_out'], header=None, names=['chrCHR', 'POS', 'END'], nrows=None)
        new['CHR'] = new['chrCHR'].str[3:]
        new['POS'] = new['POS'].astype(int)
        failed['POS'] = failed['POS'].astype(int)
        new = new.set_index(['CHR', 'POS'])
        failed = failed.set_index(['CHR', 'POS'])
        drop_coords = failed.index.intersection(ss.index)

        if not params.keep_failed_rows:
            print(len(drop_coords), 'liftover coordinates failed, will be dropped')
            ss = ss[~ss.index.isin(drop_coords)]
            ss.index = new.index
            ss.index.names = [chrCol, posCol]
            ss = ss.reset_index()[original_column_order]
        else:
            print(len(drop_coords), 'liftover coordinates failed, will be set to NA')
            ss.index.names = [chrCol, posCol]

            old_index = ss.index
            ss = ss.reset_index()[original_column_order]

            good_rows = ss.index[~old_index.isin(drop_coords)]
            print(len(good_rows))
            bad_rows = ss.index[old_index.isin(drop_coords)]
            print(len(bad_rows))

            ss.loc[good_rows, posCol] = new.index.get_level_values(1)
            print(ss)
            ss.loc[good_rows, chrCol] = new.index.get_level_values(0)
            ss.loc[bad_rows, posCol] = 0

        print(ss)
        ss.to_csv(str(output), sep='\t', header=params.has_header, index=False)
        print('End')

rule annotate_compile_suggestive_plink2:
    output:
        'Summary/plink2_all_suggestive.csv'
    input:
        expand('{cohort_dir}/Suggestive/{pheno}.txt',cohort_dir=config['plink2_cohorts'], pheno=config['plink_bin_phenos']),
        expand('{cohort_dir}/Suggestive/{pheno}.txt',cohort_dir=config['plink2_cohorts'], pheno=config['plink_quant_phenos']),
        annot='Annotations/plink2_suggestive_biofilter_genes_rsids.csv'
    resources: mem_mb=12000
    run:
        import pandas as pd
        import os

        dfs = []
        phenos = [p for p in config['plink_bin_phenos']]
        phenos.extend(config['plink_quant_phenos'])
        print(phenos)
        print(config['plink2_cohorts'])

        for c in config['plink2_cohorts']:
            for p in phenos:
                f = c + '/Suggestive/' + p + '.txt'
                if os.path.getsize(f) == 0:
                    continue
                df = pd.read_table(f, sep=' ', index_col='ID')
                df['COHORT'] = c
                df['PHENO'] = p
                dfs.append(df)

        full = pd.concat(dfs).sort_values(by=['#CHROM', 'POS'])
        print(full)
        bfDF = pd.read_csv(str(input.annot), index_col='Var_ID')
        print(bfDF)
        merged = full.merge(bfDF[['Gene', 'RSID']], left_index=True, right_index=True)
        merged.index.name = 'Var_ID'
        merged.to_csv(str(output))

rule annotate_compile_suggestive_lifted2:
    output:
        'Summary/lifted_plink2_all_suggestive.csv'
    input:
        expand('{cohort_dir}/Suggestive_Lifted/{pheno}.txt',cohort_dir=config['plink2_cohorts'], pheno=config['plink_bin_phenos']),
        expand('{cohort_dir}/Suggestive_Lifted/{pheno}.txt',cohort_dir=config['plink2_cohorts'], pheno=config['plink_quant_phenos']),
        annot='Annotations/lifted_plink2_suggestive_biofilter_genes_rsids.csv'
    resources: mem_mb=12000
    run:
        import pandas as pd
        import os

        dfs = []
        phenos = [p for p in config['plink_bin_phenos']]
        phenos.extend(config['plink_quant_phenos'])
        print(phenos)
        print(config['plink2_cohorts'])

        for c in config['plink2_cohorts']:
            for p in phenos:
                f = c + '/Suggestive_Lifted/' + p + '.txt'
                if os.path.getsize(f) == 0:
                    continue
                df = pd.read_table(f, sep=' ', index_col='ID')
                df['COHORT'] = c
                df['PHENO'] = p
                dfs.append(df)

        full = pd.concat(dfs).sort_values(by=['#CHROM', 'POS'])
        print(full)
        bfDF = pd.read_csv(str(input.annot), index_col='Var_ID')
        print(bfDF)
        merged = full.merge(bfDF[['Gene', 'RSID']], left_index=True, right_index=True)
        merged.index.name = 'Var_ID'
        merged.to_csv(str(output))

rule make_meta_input_plink2:
    output:
        'Meta/Input/{cohort_dir}.{chr}.{pheno}.sumstats.gz'
    input:
        sumstats='{cohort_dir}/Sumstats/{pheno}.sumstats.gz'
    resources: mem_mb=12000
    run:
        import pandas as pd
        import numpy as np

        df = pd.read_table(input.sumstats, index_col='ID')
        df = df[df['#CHROM'].astype(str) == str(wildcards.chr)]
        print(df)

        df = df.reset_index()

        df = df.rename(columns={'ID': 'SNP',
                                '#CHROM': 'CHR',
                                'POS': 'BP',
                                'REF': 'A2'})

        if 'BETA' not in df.columns:
            df['BETA'] = np.log(df['OR'])
            df['SE'] = (np.log(df['U95']) - np.log(df['OR'])) / 1.96

        keep_cols = ['SNP', 'BETA', 'SE', 'P', 'CHR', 'BP', 'A1', 'A2']
        df = df.drop_duplicates(subset=['CHR', 'BP', 'A2'],keep=False)
        df = df.drop_duplicates(subset=['SNP'],keep=False)
        df[keep_cols].to_csv(str(output), sep=' ', index=False)

rule make_meta_input_lifted2:
    output:
        'Meta/Input/{cohort_dir}_lifted.{chr}.{pheno}.sumstats.gz'
    input:
        sumstats='{cohort_dir}/Sumstats/{pheno}.liftover.tsv.gz'
    resources: mem_mb=12000
    run:
        import pandas as pd
        import numpy as np

        df = pd.read_table(input.sumstats, nrows=None)
        df = df[df['CHR'].astype(str) == str(wildcards.chr)]
        print(df)

        df[['chrCHR', 'pos37', 'A2', 'A1-copy']] = df['SNP'].str.split(':', expand=True)
        print(df)
        df['SNP'] = df[['chrCHR', 'BP', 'A2', 'A1']].apply(lambda x: ':'.join(x.astype(str)), axis=1)

        df['BETA'] = np.log(df['OR'])
        df['SE'] = (np.log(df['U95']) - np.log(df['OR'])) / 1.96

        keep_cols = ['SNP', 'BETA', 'SE', 'P', 'CHR', 'BP', 'A1', 'A2']
        df = df.drop_duplicates(subset=['CHR', 'BP', 'A2'], keep=False)
        df = df.drop_duplicates(subset=['SNP'], keep=False)
        df[keep_cols].to_csv(str(output), sep=' ', index=False)

rule filter_plink2_sumstats_for_rollup:
    output:
        temp('{cohort_dir}/FSS/{pheno}_{pval}.txt')
    input:
        '{cohort_dir}/Sumstats/{pheno}.sumstats.gz'
    params:
        p_col=lambda wildcards: 15 if wildcards.pheno in config['plink_quant_phenos'] else 16
    shell:
        """zcat {input} | awk '{{if (${params.p_col} <= {wildcards.pval} || NR == 1) {{print $0}}}}' > {output}"""

rule plink2_sumstats_rollup:
    output:
        '{cohort_dir}/Sumstats/all_phenos_p_cutoff_{pval}.csv.gz'
    input:
        quant_sumstats = expand('{{cohort_dir}}/FSS/{pheno}_{{pval}}.txt', pheno=config['plink_quant_phenos']),
        bin_sumstats = expand('{{cohort_dir}}/FSS/{pheno}_{{pval}}.txt', pheno=config['plink_bin_phenos'])
    resources: mem_mb=12000
    run:
        import pandas as pd

        dfs = []
        phenos = config['plink_quant_phenos']
        phenos.extend(config['plink_bin_phenos'])
        pval = wildcards.pval
        c = wildcards.cohort_dir
        for p in phenos:
            f = c + '/FSS/' + p + '_' + pval + '.txt'
            if os.path.getsize(f) == 0:
                continue
            df = pd.read_table(f, sep='\t')
            print(df)
            df['COHORT'] = c
            df['PHENO'] = p
            dfs.append(df)

        full = pd.concat(dfs).sort_values(by=['#CHROM', 'POS'])
        full.to_csv(str(output), index=False)

rule plink2_cohort_sumstats_rollup:
    output:
        'Summary/all_plink2_cohorts_all_phenos_p_cutoff_{pval}.csv.gz'
    input:
        expand('{cohort_dir}/Sumstats/all_phenos_p_cutoff_{{pval}}.csv.gz', cohort_dir=config['plink2_cohorts'])
    resources: mem_mb=30000
    run:
        import pandas as pd

        dfs = []
        for f in input:
            if os.path.getsize(f) == 0:
                continue
            df = pd.read_csv(f)
            dfs.append(df)

        full = pd.concat(dfs).sort_values(by=['#CHROM', 'POS'])
        full.to_csv(str(output),index=False)
