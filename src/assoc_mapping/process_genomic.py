"""
This script process the "genomic" output from populations,
turning the genotype encoding into a proportion of homozygosity
and removing SNPs that are homozygous or missing in mothers
from their respective offspring. Since genomics output file from
populations are typically huge, the script takes quite a long time
to run although it is parallelized and will automatically use all
core available on the host. It takes 2 arguments: the path to the
folder containing populations files, used as input, and the folder where
the output will be written.
"""
# Cyril Matthey-Doret
# 11.08.2017

from os import path, walk
import re
import argparse
import numpy as np
import pandas as pd
from multiprocessing import Pool  # Parallel computing support
from functools import partial  # "freeze" arguments when mapping function


# ========== PARSING COMMAND LINE ARGUMENTS ==========#

parser = argparse.ArgumentParser(description="This script processes the \
                                 'genomic' output from STACKS's populations \
                                 module to compute the proportion of \
                                 homozygous individuals at each genomic \
                                 position.")
parser.add_argument('pop_files', type=str,
                    help='Path to the folder containing \
                                 populations output files. Used as input')
parser.add_argument('out', type=str,
                    help='Folder where output will be written')
parser.add_argument('--keep_all', action='store_true',
                    help='Keep all SNPs, even if missing or\
                                 homozygous in the mother.')
parser.add_argument('--pool_output', action='store_true',
                    help='Pool output instead of computing \
                                 proportions per family')
parser.add_argument('--sample_list', type=str,
                    help='File containing the list of samples with family, \
                    name, sex and generation information. Defaults to \
                    "data/individuals.tsv".', default='data/individuals.tsv')

args = parser.parse_args()

# ========== DEFINING FUNCTIONS ==========#


def unify_genomic(pop_path, pop):
    """
    Reads populations "genomic" output files from each family's folder and
    gathers them into a single pandas DataFrame object.
    :param path: Path containing the family subfolders containing populations
    output files.
    :param pop: dataframe containing name of all individuals in correct final
    order.
    :returns: A DataFrame containing all sites in all individuals.
    """

    # pop['real_idx'] = pop.index
    # Adding trailing slash if not provided
    if pop_path[-1] != '/':
        pop_path += '/'
    first_fam = True
    # Iterating over family subfolders
    for subdir, dirs, files in walk(pop_path):
        if path.basename(subdir):
            # Read file from subfolders
            sum_file = path.join(subdir, "batch_1.sumstats.tsv")
            with open(sum_file, 'r') as f:
                # Count population rows starting at -1 to exclude header line
                pops = -1
                for line in f:
                    if re.search('^#', line):
                        pops += 1
            # Reading only population rows
            fam_samples = pd.read_csv(sum_file, sep='\t',
                                      nrows=pops, header=None)
            # Get individual names in original order
            f_names_nest = [fam_samples.iloc[:, 1][n].split(',') for n in
                            range(pops)]
            # Flattening nested lists while maintaining order
            fam_names = [item for nest in f_names_nest for item in nest]
            fam_names = pd.DataFrame({'Name': fam_names})
            # Get final index of names
            tmp_pop = fam_names.merge(pop, on='Name', how='inner')
            # Incrementing index by 3 since there is 3 columns before indv
            # tmp_idx = tmp_idx.real_idx + 3
            tmp = pd.read_csv(path.join(subdir, "batch_1.genomic.tsv"),
                              sep='\t', header=None, skiprows=1)
            # Rename sample columns with final indices
            tmp.rename(columns={x: tmp_pop.Name[x-3] for x in
                                range(3, tmp.shape[1])}, inplace=True)
            if first_fam:
                # Assign dataframe on first iteration only
                united = tmp
                first_fam = False
            else:
                # Graft each family's samples onto final df as new columns
                # Using outer merge on chromosome, bp and Locus ID
                united = united.merge(tmp, how='outer', on=range(3))
    united = united.fillna(0)  # Change all NAs to the "missing" code
    return united


def gen_decode(encoded):
    """"
    This function decodes numeric genotypes and
    replaces it with E (heterozygous), O (homozygous)
    or M (missing).
    :returns: A dataframe containing the genotype letters.
    """

    genodict = {}
    for code in range(11):
        # Genotype encoding in genomic output from populations:
        # Missing bases are encoded as 0
        # Homozygous genotypes are: 1,5,8,10
        # Heterozygous genotypes are all others (except 0)
        # Building dictionary for numeric genotype translation
        if code in [1, 5, 8, 10]:  # Homozygous codes
            genodict[code] = ''.join(['O',str(code)])
        elif code == 0:  # Missing
            genodict[code] = 'M'
        else:
            genodict[code] = 'E'  # All others are heterozygous
    # Rarely, rows are filled this value. I assume this is a STACKS issue.
    decoded = encoded.applymap(lambda r: genodict.get(r, 'M'))
    return decoded


def mother_hom(geno, pop):
    """
    This function runs on a numpy array that has already been
    transformed with gen_decode and sets sites that are homozygous/missing in
    mothers to missing in their whole family. If the mother is not available,
    sites that are homozygous or missing in all offspring in the family are
    used instead as a proxy.
    :param pop: a dataframe containing individual names and their respective
    families. The names need to be in the same order as the columns in geno.
    :param geno: a numpy array that will be processed
    """

    for f in np.unique(pop.Family):  # Iterate over mothers
        fam = pop.loc[pop.Family == f, :]  # Subssetting samples from family
        # Get mother idx
        mother_name = fam.Name[fam.Generation == 'F3'].tolist()
        # hom/missing mother sites
        fam_SNP = np.where(geno.loc[:, mother_name] != 'E')[0]
        if not mother_name:  # If the mother is not available
            # Use sites where no individual in the fam is heterozygous instead
            fam_df = pd.DataFrame(geno[fam.Name.tolist()].T)
            # Count number of diff genotypes among offspring (E, M, O1, O2)
            fam_SNP = fam_df.apply(lambda x: len(x.unique()))
            # If at least one offspring is het: SNP used automatically
            fam_SNP += (fam_df == "E").sum() * 3
            # If missing among genotypes, decrease number by one for the SNP
            fam_SNP -= (fam_df == 'M').any()

            # Include SNP if number of genotypes > 1 (het = auto include)
            fam_SNP = fam_SNP[fam_SNP < 2].index.tolist()

        # Change those sites to M in all indv with same family
        geno.loc[fam_SNP, fam.Name] = 'M'
        # If using pandas <0.19, the above line will fail. The loop below can
        # be used instead, but is much slower
        # for snp in fam_SNP:
        #    for family in fam.Name:
        #        geno.set_value(snp, family, 'M')
    return geno


def prop_hom(pop, geno):
    """
    This function computes the proportion of homozygous individuals (cols) at
    each SNP (row) in a numpy array containing decoded allelic state (O,E,M).
    It computes this proportion both by sex, and in all individuals.
    :param geno: Pandas DataFrame with sites as rows and individuals as cols.
    :param pop: Dataframe containing the sex of each individual and its name.
    :returns: a Pandas DataFrame object with the proportion of homozygous
    females, males and all individuals at each site and the number of
    individuals where it was present.
    """
    # Suppresses warning when numpy divides by 0
    np.seterr(divide='ignore', invalid='ignore')
    # Number of males and females
    # N = {sex:pop.Sex[pop.Sex == sex].shape[0] for sex in ['M','F']}
    N = {'M': pop.Sex[pop.Sex == 'M'].shape[0],
         'F': pop.Sex[pop.Sex == 'F'].shape[0]}
    # Get sample names by sex
    sex_id = {'M': pop.Name[pop.Sex == 'M'],
              'F': pop.Name[pop.Sex == 'F']}
              
    # Un-specify the genotype of homozygous individuals (to comptabilise recomb)
    geno = geno.replace(['O1','O5','O8', 'O10'], "O")
    # Counting how many individuals are used to compute proportion at each SNP
    sample_size = {}  # Number of individuals in which each site is found
    hom = {}  # proportion of individuals in which each site is homozygous
    for sex in N:
        # Looping over sexes
        dff = {}
        for t in ['O', 'M', 'E']:
            dff[t] = (geno.loc[:, sex_id[sex].values] == t).T.sum().astype(float)
        sample_size[sex] = dff['E']+dff['O']
        hom[sex] = np.divide(dff['O'], (dff['O'] + dff['E']))

    # Building output dataframe with all relevant stats
    out_df = pd.DataFrame({
        "N.Samples": sample_size['F'] + sample_size['M'],
        "Prop.Hom": ((sample_size['M'] * hom['M'].fillna(0) +
                      sample_size['F'] * hom['F'].fillna(0)) /
                     (sample_size['F'] + sample_size['M'])).round(3),
        "N.Males": sample_size['M'],
        "N.Females": sample_size['F'],
        "Prop.Hom.F": hom['F'].round(3),
        "Prop.Hom.M": hom['M'].round(3)
        })
    return out_df


def split_fam_prop(df, pop, parallel=True):
    """
    This function is intended as a wrapper for prop_hom, so that it will only
    compute stats per family and return family info associated with every stat.
    :param df: Pandas DataFrame with sites as rows and individuals as columns.
    :param pop: Dataframe containing the sex, name and Family of each
    individual in the same order as the df individuals columns.
    :param parallel: Boolean value. Should the script exploit multiple cores if
    available ?
    :returns: a Pandas DataFrame object with an added family column and at each
    site for each family, the proportion of homozygous females, males and all
    individuals and the number of individuals where it was present.
    """

    df_list = []  # List to contain each family's df
    for fam in pop.Family.unique():
        # Iterating over families and subsetting df for each
        fam_id = pop.Name[pop.Family == fam]
        genofam = df.loc[:, fam_id]
        popfam = pop.loc[pop.Name.isin(fam_id), :]
        if parallel:
            fam_df = parallel_func(prop_hom, genofam, f_args=(popfam,))
        else:
            fam_df = prop_hom(popfam, genofam)
        fam_df["Family"] = fam
        df_list.append(fam_df)

    return pd.concat(df_list)


def parallel_func(f, df, f_args=[], chunk_size=1000):
    """
    Parallelizes a function that runs on a dataframe by splitting the dataframe
    into small chunks by rows and distributing chunks across several processes.
    :param f: Target function that will be parallelized
    :param df: pandas dataframe to be used as input
    :param f_args: optional arguments for the function to be parallelized. Need
    to be an iterable (list or tuple).
    :param chunk_size: size of the chunks in which df is split. Default=1000
    :returns: the processed dataframe reconstructed by combining output from
    all processes
    """

    # Create pool of processes, size depends on number of core available
    pool = Pool(processes=4)
    tot_rows = df.shape[0]
    chunks = range(0, tot_rows, chunk_size)  # Start positions of chunks
    # Split df into chunks
    chunked_df = [df.iloc[c: (c+min(chunk_size, tot_rows)), ] for c in chunks]
    func = partial(f, *f_args)  # Unpacking optional fixed arguments.
    result = pool.map(func, chunked_df)  # Mapping function to chunks.
    # Concatenating into single df. Order is preserved
    pool.terminate()

    return pd.concat(result)

# ========== LOADING AND PROCESSING DATA ==========#
# Path to STACKS populations folder and output file


in_path = args.pop_files
out_prefix = 'grouped_'
if args.pool_output:
    out_prefix += "outpool_"
if args.keep_all:
    out_prefix += "keep_"
out_path = path.join(args.out, (out_prefix + "prophom.tsv"))
indv_path = args.sample_list  # family and sex information
indv = pd.read_csv(indv_path, sep='\t')  # Family and sex info
# Preparing data structure to match sample names and families with columns


try:
    # Names in correct order
    samples = pd.read_csv(path.join(in_path, "batch_1.sumstats.tsv"),
                          sep='\t', nrows=2, header=None)

    # Concatenating populations
    names = samples.iloc[:, 1][0].split(',') + \
        samples.iloc[:, 1][1].split(',')
except pd.errors.ParserError:
    # In case only 1 sex is present
    samples = pd.read_csv(path.join(in_path, "batch_1.sumstats.tsv"),
                          sep='\t', nrows=1, header=None)

    # Concatenating populations
    names = samples.iloc[:, 1][0].split(',')

names = pd.DataFrame({'Name': names})
# Adding family and sex, keeping order
pop = names.merge(indv, on='Name', how='left')
genomic = pd.read_csv(path.join(in_path, "batch_1.genomic.tsv"),
                      sep='\t', header=None, skiprows=1)
# select only samples cols and reindex from 0
gen_indv = genomic.iloc[:, 3:].T.reset_index(drop=True).T
# Replacing numeric header with corresponding sample name
gen_indv.rename(columns=lambda x: pop.Name[x], inplace=True)
print("Processing {0} sites across {1} samples.".format(gen_indv.shape[0],
                                                        gen_indv.shape[1]))
print("files loaded")

# ========== RUNNING CODE ==========#
# Decoding numeric genotypes into states (het, hom, missing)
state = parallel_func(gen_decode, gen_indv)
print("genotypes decoded")
# Will run unless user explicitly set the --keep_all parameter
if not args.keep_all:
    # Remove SNPs that are hom./missing in mothers from their family
    state = mother_hom(state, pop)
    print("Mother homozygous and missing sites removed")
# Computing proportion of homozygous indv at each site
if args.pool_output:
    prop = prop_hom(pop, state)
else:
    prop = split_fam_prop(state, pop, parallel=False)
print("homozygosity stats calculated")

# ========== SAVING OUTPUT ==========#
# Merging Chromosomal positions with proportion of homozygosity into 1 df
prop = genomic.iloc[:, 0: 3].merge(prop, left_index=True, right_index=True)
state = genomic.iloc[:, 0: 3].merge(state, left_index=True, right_index=True)
prop.rename(columns={0: "Locus.ID", 1: "Chr", 2: "BP"}, inplace=True)
state.rename(columns={0: "Locus.ID", 1: "Chr", 2: "BP"}, inplace=True)
state_path = path.join(args.out,
                       (out_prefix.replace("outpool_", "") + "geno_EOM.tsv"))
state = state.replace(['O1','O5','O8', 'O10'], "O")
state.to_csv(state_path, sep='\t', index=False, na_rep='NA')
prop.to_csv(out_path, sep='\t', index=False, na_rep='NA')
print("Output saved to {0}".format(out_path))
