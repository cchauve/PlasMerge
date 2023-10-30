#!/usr/bin/env python

'''
Importing data from NCBI
'''

import sys
import os
import argparse
import pandas as pd
import xml.etree.ElementTree as ET
import subprocess
from collections import defaultdict

## Dataframes functions

# Separators
NCBI_SEP = '; '
LOG_SEP = ','
PD_SEP = ':'
ENTREZ_SEP = ','
# Column names of dataframe
BIOSAMPLE_KEY = 'Assembly BioSample Accession'
ASSEMBLY_KEY = 'Assembly Accession'
SRA_KEY = 'SRA Accession'
GB_CHR_KEY = 'GenBank Chromosome Accession'
GB_PLS_KEY = 'GenBank Plasmid Accession'
GB_LEN_KEY = 'Molecule Length'
# Default directory where to write BioSamples XML files
BIOSAMPLES_OUT_DIR = 'BioSamples'
# Default directory where to write GenBank docsum files
GENBANK_OUT_DIR = 'GenBank'
# Default directory where to write SRA docsum files
SRA_OUT_DIR = 'SRA'


def read_tsv_file(in_file):
    '''
    Read a TSV file and returns a dataframe
    '''
    return pd.read_table(in_file)

def df_extract_column(in_df, in_col):
    '''
    Extract from a dataframe a column of name in_col
    Returns dict(index -> value)
    '''
    out_dict = {}
    for idx,row in in_df.iterrows():
        value = row[in_col]
        out_dict[idx] = value        
    return out_dict

def df_add_column(in_df, in_dict, in_col, join=False):
    '''
    Add a column to in_df from dictionary in_dict indexed by in_df index
    '''
    for idx,row in in_df.iterrows():
        if len(in_dict[idx])>0:
            if join:
                in_df.at[idx, in_col] = NCBI_SEP.join(in_dict[idx])
            else:
                in_df.at[idx, in_col] = in_dict[idx]
        else:
            in_df.at[idx, in_col] = pd.NA

def df_expand_column(in_df, in_dict, in_col, join=False):
    '''
    Expand an existing column to in_df from dictionary in_dict indexed by in_df index
    Assumption:
    if len(in_dict[idx])>0 the in_df.at[idx, in_col] is not NaN
    '''
    for idx,row in in_df.iterrows():
        if len(in_dict[idx]) > 0:
            expanded_list = in_df.at[idx, in_col].split(NCBI_SEP) + in_dict[idx]
            if join:
                in_df.at[idx, in_col] = NCBI_SEP.join(expanded_list)
            else:
                in_df.at[idx, in_col] = expanded_list

def df_export_column(in_df, in_col, out_file_path, sep=',', filter_fun=lambda x: x):
    '''
    Export the values of a column to a TSV file indexed by biosample
    '''
    with open(out_file_path, 'w') as out_file:
        for idx,row in in_df.iterrows():
            biosample = row[BIOSAMPLE_KEY]
            value = filter_fun(row[in_col])
            if value is not None:
                out_file.write(f'{biosample}\t{sep.join(value.split(NCBI_SEP))}\n')

def df_check_duplicates(in_df, in_col):
    '''
    Returns a list of values that appears more than once in the column in_col
    '''
    result = []
    in_df_grouped = in_df.groupby([in_col])
    for group_idx,group_df in in_df_grouped:
        if len(group_df.index) > 1:
            result.append(group_idx)
    return result

def df_delete_rows(in_df, in_col, values_list):
    '''
    Delete (in_place) rows whose in_col value is in values_list
    '''
    in_df.drop(in_df[in_df[in_col].isin(values_list)].index, inplace=True)
    
                
## Docsum format functions

DOCSUM_ERROR_KEY = 'DOCSUM_ERROR'
DOCSUM_MISSING_FILE_KEY = 'DOCSUM_MISSING'
DOCSUM_BIOSAMPLE_ID_ERROR_KEY = 'ID_ERROR'
DOCSUM_BIOSAMPLE_ID_MISSING_KEY = 'ID_MISSING'
DOCSUM_BIOSAMPLE_ID_AMBIGUOUS_KEY = 'ID_AMBIGUOUS'

def docsum_get_acc_text(docsum_file, acc_key):
    '''
    Read a docsum file as a text file.
    Returns a list of the accession IDs for key acc_key
    Used for incorrectly formetted docsum files
    '''
    accessions_ids = []
    with open(docsum_file) as in_file:
        for docsum_line in in_file.readlines():
            line = docsum_line.strip()
            if line.startswith(f'<{acc_key}'):
                acc_id = line.split('>')[1].split('<')[0]
                if acc_id not in accessions_ids:
                    accessions_ids.append(acc_id)
    return accessions_ids

def docsum_get_acc_xml(docsum_file, acc_key):
    '''
    Read a docsum file.
    Returns a list of the accssion IDs for database acc_key
    If the file is not in correct docsum format, return None
    '''
    try:
        accessions_ids = []
        docsum_tree = ET.parse(docsum_file)
        entries = docsum_tree.getroot().findall('DocumentSummary')
        for entry in entries:
            for acc_tag in entry.iter(acc_key):
                acc_id = acc_tag.text
                if acc_id not in accessions_ids:
                    accessions_ids.append(acc_id)
        return accessions_ids
    except ET.ParseError:
        print(f'ERROR\treading {docsum_file}', file=sys.stderr)
        return None

def check_biosample_docsum_files(in_df, out_dir, key):
    '''
    Returns a dictionary of encountered errors for each BioSample in in_df
    - DOCSUM_MISSING_FILE_KEY: missing DOCSUM file
    - DOCSUM_ERROR_KEY: incorrectly formatted DOCSUM file
    - DOCSUM_BIOSAMPLE_ID_ERROR_KEY: BioSample ID in file different from input BioSample ID
    - DOCSUM_BIOSAMPLE_ID_MISSING_KEY: No BioSample ID in file
    - DOCSUM_BIOSAMPLE_ID_AMBIGUOUS_KEY: more than one BioSample ID in file 
    '''
    errors = defaultdict(list)
    for idx,row in in_df.iterrows():
        biosample = row[BIOSAMPLE_KEY]
        docsum_file = os.path.join(out_dir, f'{biosample}.docsum')
        if os.path.isfile(docsum_file):
            docsum_error = False
            biosample_check = docsum_get_acc_xml(docsum_file, key)
            if biosample_check is None:
                biosample_check = docsum_get_acc_text(docsum_file, key)
                docsum_error = True
            if docsum_error:
                errors[idx].append(DOCSUM_ERROR_KEY)
            if len(biosample_check) > 1:
                errors[idx].append(DOCSUM_BIOSAMPLE_ID_AMBIG_KEY)
            elif len(biosample_check) == 0:
                errors[idx].append(DOCSUM_BIOSAMPLE_ID_MISSING_KEY)
            elif biosample_check[0] != biosample:
                errors[idx].append(DOCSUM_BIOSAMPLE_ID_ERROR_KEY)
        else:
            errors[idx].append(DOCSUM_MISSING_FILE_KEY)
    return errors

## Intermediate, BioSample-indexed, files

def create_biosample_docsum_file(biosample, docsum_file, max_attempts=5):
    docsum_cmd = ['efetch', '-db', 'biosample', '-id', biosample, '-format', 'docsum']
    with open(docsum_file, 'w') as out_file:
        for attempt in range(max_attempts):
            try:
                subprocess.run(docsum_cmd, check=True, stdout=out_file)
                success = True
                break
            except subprocess.CalledProcessError:
                success = False
    if not success:
        print(f'ERROR\tcreating {docsum_file}', file=sys.stderr)

def create_biosample_docsum_files(in_df, out_dir):
    '''
    Create, for BioSample accessions IDs in in_df, one BioSample docsum file per ID.
    File names: out_dir/ID.docsum
    Notes:
    We could try to use the Entrez python API available in BioPython but it seems to not read properly BioSample IDs.
    '''
    os.makedirs(out_dir, exist_ok=True)
    errors = defaultdict(list)
    for idx,row in in_df.iterrows():
        biosample = row[BIOSAMPLE_KEY]
        docsum_file = os.path.join(out_dir, f'{biosample}.docsum')
        create_biosample_docsum_file(biosample, docsum_file)

def get_tag_biosample_docsum_files(in_df, in_dir, tag):
    '''
    Extract from a docsum file the value of a tag defied by its prefix
    Returns a directory index -> list of tag values
    '''
    out_dict = defaultdict(list)
    for idx,row in in_df.iterrows():
        biosample = row[BIOSAMPLE_KEY]
        docsum_file = os.path.join(in_dir, f'{biosample}.docsum')
        tag_values = docsum_get_acc_text(docsum_file, tag)
        out_dict[idx] = tag_values
    return out_dict
        
def get_SRA_biosample_docsum_files(in_df, in_dir):
    '''
    Read BioSample docsum files for the biosamples in in_df, located in 
    directory in_dir and extracts for each the associated SRA.
    Returns a dictionary index -> list of SRA IDs
    '''
    return get_tag_biosample_docsum_files(in_df, in_dir, 'Id db="SRA"')

def get_Organism_biosample_docsum_files(in_df, in_dir):
    '''
    Read BioSample docsum files for the biosamples in in_df, located in 
    directory in_dir and extracts for each the associated Organism.
    Returns a dictionary index -> [Organism]
    '''
    return get_tag_biosample_docsum_files(in_df, in_dir, 'OrganismName>')


def create_genbank_docsum_file(assembly, docsum_file, max_attempts=5):
    '''
    Creates a GenBank file from an assembly ID
    Attempts at most max_attempts to create the file
    '''
    docsum_cmd_esearch = ['esearch', '-db', 'assembly', '-query', assembly]
    docsum_cmd_elink = ['elink', '-target', 'nucleotide', '-name', 'assembly_nuccore_insdc']
    docsum_cmd_efetch = ['efetch', '-format', 'docsum']
    with open(docsum_file, 'w') as out_file:
        for attempt in range(max_attempts):
            try:
                esearch_process = subprocess.run(docsum_cmd_esearch, check=True, capture_output=True)
                elink_process = subprocess.run(docsum_cmd_elink, input=esearch_process.stdout, check=True, capture_output=True)
                efetch_process = subprocess.run(docsum_cmd_efetch, input=elink_process.stdout, check=True, stdout=out_file)
                success = True
                break
            except subprocess.CalledProcessError:
                success = False
    if not success:
        print(f'ERROR\tcreating {docsum_file}', file=sys.stderr)

            
def create_genbank_docsum_files(in_df, out_dir):
    '''
    Create, one GenBank docsum file per BioSample ID.
    File names: out_dir/ID.docsum
    '''
    os.makedirs(out_dir, exist_ok=True)
    for idx,row in in_df.iterrows():
        assembly,biosample = row[ASSEMBLY_KEY],row[BIOSAMPLE_KEY]
        docsum_file = os.path.join(out_dir, f'{biosample}.docsum')
        create_genbank_docsum_file(assembly, docsum_file)

def get_GB_genbank_docsum_files(in_df, in_dir):
    '''
    Read GenBank docsum files for the biosamples in in_df, located in 
    directory in_dir and extracts for each the GenBank accessions ID
    for chromosomes and plasmids and their lengths
    Returns a dictionary index -> dict{'chromosome': list ID, 'plasmid': list ID, 'len': list str(ID:length)}
    Assumption:
    Any molecule without plasmid in them is a chromosome
    '''
    out_dict = {}
    for idx,row in in_df.iterrows():
        biosample = row[BIOSAMPLE_KEY]
        docsum_file = os.path.join(in_dir, f'{biosample}.docsum')
        out_dict[idx] = {'chromosome': [], 'plasmid': [], 'len': []}
        try:
            xml_tree = ET.parse(docsum_file)
            xml_entries = xml_tree.getroot().findall('DocumentSummary')
            for entry in xml_entries:
                title = entry.find('Title').text.lower()
                if 'plasmid' in title: molecule = 'plasmid'
                else: molecule = 'chromosome'
                extra = entry.find('Extra').text.split('|')
                if 'emb' in extra: acc_id = extra[extra.index('emb')+1]
                elif 'gb' in extra: acc_id = extra[extra.index('gb')+1]
                elif 'dbj' in extra: acc_id = extra[extra.index('dbj')+1]
                else: acc_id = None
                length = int(entry.find('Slen').text)
                if acc_id is not None and acc_id not in out_dict[idx][molecule]:
                    out_dict[idx][molecule].append(acc_id)
                    out_dict[idx]['len'].append(f'{acc_id}{PD_SEP}{length}')
        except ET.ParseError:
            print(f'ERROR\treading {docsum_file}', file=sys.stderr)
    return out_dict

def create_SRA_docsum_file(SRA, docsum_file, max_attempts=5):
    '''
    Creates a GenBank file from an assembly ID
    Attempts at most max_attempts to create the file
    '''
    docsum_cmd = ['efetch', '-db', 'sra', '-id', SRA, '-format', 'docsum']
    with open(docsum_file, 'w') as out_file:
        for attempt in range(max_attempts):
            try:
                subprocess.run(docsum_cmd, check=True, stdout=out_file)
                success = True
                break
            except subprocess.CalledProcessError:
                success = False
    if not success:
        print(f'ERROR\tcreating {docsum_file}', file=sys.stderr)
            
def create_SRA_docsum_files(in_df, out_dir):
    '''
    Create, one SRA docsum file per SRA ID.
    File names: out_dir/ID.docsum
    '''
    os.makedirs(out_dir, exist_ok=True)
    for idx,row in in_df.iterrows():
        SRA,biosample = row[SRA_KEY],row[BIOSAMPLE_KEY]
        if pd.notnull(SRA):
            docsum_file = os.path.join(out_dir, f'{biosample}.docsum')
            create_SRA_docsum_file(SRA, docsum_file)

def get_run_SRA_docsum_files(in_df, in_dir):
    '''
    Read SRA docsum files for the biosamples in in_df, located in 
    directory in_dir and extracts for each the list of SRA runs
    Returns a dictionary index -> list (instrument:run acc)
    '''
    out_dict = {}
    for idx,row in in_df.iterrows():
        biosample = row[BIOSAMPLE_KEY]
        docsum_file = os.path.join(in_dir, f'{biosample}.docsum')
        out_dict[idx] = []
        if os.path.isfile(docsum_file):
            try:
                xml_tree = ET.parse(docsum_file)
                xml_entries = xml_tree.getroot().findall('DocumentSummary')
                for entry in xml_entries:
                    platform = entry.find('ExpXml').find('Summary').find('Platform').text.upper()
                    runs = entry.find('Runs').findall('Run')
                    for run in runs:
                        run_acc = run.attrib['acc']
                        out_dict[idx].append(f'{platform}{PD_SEP}{run_acc}')
            except ET.ParseError:
                print(f'ERROR\treading {docsum_file}', file=sys.stderr)
    return out_dict

def entrez_download_file(in_db, in_query, in_fmt, out_file_path, max_attempts=5):
    '''
    Download the data for accessions in_query, from database in_db in format in_fmt
    into file out_file_path
    '''
    esearch_cmd = ['esearch', '-db', in_db, '-query', in_query]
    efetch_cmd  = ['efetch', '-format', in_fmt]
    with open(out_file_path, 'w') as out_file:
        for attempt in range(max_attempts):
            try:
                esearch_process = subprocess.run(esearch_cmd, check=True, capture_output=True)
                efetch_process = subprocess.run(efetch_cmd, input=esearch_process.stdout, check=True, stdout=out_file)
                success = True
                break
            except subprocess.CalledProcessError:
                success = False
        if not success:
            print(f'ERROR\tdownloading {out_file_path}', file=sys.stderr)
        else:            
            print(f'SUCCESS\tdownloading {out_file_path}', file=sys.stdout)
            
## Log functions

def write_log(log_file, access, in_df, key, header, errors):
    with open(log_file, access) as out_log:
        for idx,error in errors.items():
            if len(error)>0:
                acc = in_df.at[idx, key]
                out_log.write(f'{header}\t{acc}\t{LOG_SEP.join(error)}\n')

## Commands functions

def cmd_genomes(in_tsv_file, out_tsv_file, out_log_file):
    '''
    Reads an input genome files and augment it by SRA and GenBank IDs for plasmids
    '''
    genomes_df = read_tsv_file(in_tsv_file)

    # Adding SRA information
    create_biosample_docsum_files(genomes_df, BIOSAMPLES_OUT_DIR)
    biosamples_docsum_errors = check_biosample_docsum_files(genomes_df, BIOSAMPLES_OUT_DIR, 'Accession')
    write_log(out_log_file, 'w', genomes_df, BIOSAMPLE_KEY, 'BioSample', biosamples_docsum_errors)
    Organism_dict = get_Organism_biosample_docsum_files(genomes_df, BIOSAMPLES_OUT_DIR)
    df_add_column(genomes_df, Organism_dict, 'Organism', join=True)
    SRA_dict = get_SRA_biosample_docsum_files(genomes_df, BIOSAMPLES_OUT_DIR)
    df_add_column(genomes_df, SRA_dict, 'SRA Accession', join=True)
    create_SRA_docsum_files(genomes_df, SRA_OUT_DIR)
    SRA_docsum_errors = check_biosample_docsum_files(genomes_df, SRA_OUT_DIR, 'BioSample')
    write_log(out_log_file, 'a', genomes_df, BIOSAMPLE_KEY, 'SRA', SRA_docsum_errors)
    SRA_exp_dict = get_run_SRA_docsum_files(genomes_df, SRA_OUT_DIR)
    df_expand_column(genomes_df, SRA_exp_dict, 'SRA Accession', join=True)
    # Adding GenBank information
    create_genbank_docsum_files(genomes_df, GENBANK_OUT_DIR)
    genbank_docsum_errors = check_biosample_docsum_files(genomes_df, GENBANK_OUT_DIR, 'BioSample')
    write_log(out_log_file, 'a', genomes_df, BIOSAMPLE_KEY, 'GenBank', genbank_docsum_errors)    
    GenBank_dict = get_GB_genbank_docsum_files(genomes_df, GENBANK_OUT_DIR)
    GenBank_dict_chromosome = {idx: acc_ids['chromosome'] for idx, acc_ids in GenBank_dict.items()}
    df_add_column(genomes_df, GenBank_dict_chromosome, GB_CHR_KEY, join=True)
    GenBank_dict_plasmid = {idx: acc_ids['plasmid'] for idx, acc_ids in GenBank_dict.items()}
    df_add_column(genomes_df, GenBank_dict_plasmid, GB_PLS_KEY, join=True)
    GenBank_dict_length = {idx: acc_ids['len'] for idx, acc_ids in GenBank_dict.items()}
    df_add_column(genomes_df, GenBank_dict_length, GB_LEN_KEY, join=True)
    genomes_df.to_csv(out_tsv_file, sep='\t', index=False)

def cmd_entrez_download(in_tsv_file, in_col, in_db, out_fmt, out_ext, out_dir):
    '''
    Reads an input TSV file, extract the accessions ID from column in_col.
    Fetch the corresponding data from database in_db, in format out_fmt, written in a file
    out_dir/{biosample}.{out_ext}
    '''
    in_df = read_tsv_file(in_tsv_file)
    biosamples = df_extract_column(in_df, BIOSAMPLE_KEY)
    accessions = df_extract_column(in_df, in_col)
    os.makedirs(out_dir, exist_ok=True)
    for idx,acc in accessions.items():
        if pd.notnull(acc):
            biosample = biosamples[idx]
            out_file_path = os.path.join(out_dir, f'{biosample}.{out_ext}')
            in_query = ENTREZ_SEP.join(acc.split(NCBI_SEP))
            entrez_download_file(in_db, in_query, out_fmt, out_file_path)            

def cmd_delete_entries(in_tsv_file, in_col, in_values, out_tsv_file=None):
    in_df = read_tsv_file(in_tsv_file)
    df_delete_rows(in_df, in_col, in_values)
    if out_tsv_file is None:
        in_df.to_csv(in_tsv_file, sep='\t', index=False)
    else:
        in_df.to_csv(out_tsv_file, sep='\t', index=False)

def cmd_export_SRA(in_df, platform, out_tsv_file):
    '''
    Export the SRA Accession column keeping only the accession for the
    sequencing platform given by platform.
    '''
    def filter_fun(value):
        if pd.isnull(value): return None
        SRA_list = [
            x.split(PD_SEP)[1]
            for x in value.split(NCBI_SEP)[1:]
            if x.split(PD_SEP)[0]==platform
        ]
        if len(SRA_list) == 0: return None
        else: return NCBI_SEP.join(SRA_list)
    df_export_column(
        in_df, SRA_KEY, out_tsv_file, filter_fun=filter_fun
    )

def main():
    parser = argparse.ArgumentParser(description='Importing NCBI data.')
    subparsers = parser.add_subparsers(help='sub-command help')
    # genomes command arguments
    genomes_parser = subparsers.add_parser('genomes')
    genomes_parser.set_defaults(cmd='genomes')
    genomes_parser.add_argument('input', type=str, help='Input TSV file')
    genomes_parser.add_argument('output', type=str, help='Output TSV file')
    genomes_parser.add_argument('log', type=str, help='Log file')
    # entrez_download command arguments
    entrez_download_parser = subparsers.add_parser('entrez_download')
    entrez_download_parser.set_defaults(cmd='entrez_download')
    entrez_download_parser.add_argument('input', type=str, help='Input TSV file')
    entrez_download_parser.add_argument('column', type=str, help='Input TSV file column to read')
    entrez_download_parser.add_argument('db', type=str, help='Entrez database')
    entrez_download_parser.add_argument('format', type=str, help='File format')    
    entrez_download_parser.add_argument('ext', type=str, help='File extension')    
    entrez_download_parser.add_argument('out_dir', type=str, help='Output directory')    
    # export_column command arguments
    export_column_parser = subparsers.add_parser('export_column')
    export_column_parser.set_defaults(cmd='export_column')
    export_column_parser.add_argument('input', type=str, help='Input TSV file')
    export_column_parser.add_argument('column', type=str, help='Input TSV file column to read')
    export_column_parser.add_argument('out_file', type=str, help='Output TSV file')
    # export_SRA command arguments
    export_SRA_parser = subparsers.add_parser('export_SRA')
    export_SRA_parser.set_defaults(cmd='export_SRA')
    export_SRA_parser.add_argument('input', type=str, help='Input TSV file')
    export_SRA_parser.add_argument('platform', type=str, help='Sequencing platform')
    export_SRA_parser.add_argument('out_file', type=str, help='Output TSV file')    
    # check_duplicates command arguments
    check_duplicates_column_parser = subparsers.add_parser('check_duplicates')
    check_duplicates_column_parser.set_defaults(cmd='check_duplicates')
    check_duplicates_column_parser.add_argument('input', type=str, help='Input TSV file')
    check_duplicates_column_parser.add_argument('column', type=str, help='Input TSV file column to read')
    # delete_entries command arguments
    delete_entries_column_parser = subparsers.add_parser('delete_entries')
    delete_entries_column_parser.set_defaults(cmd='delete_entries')
    delete_entries_column_parser.add_argument('input', type=str, help='Input TSV file')
    delete_entries_column_parser.add_argument('column', type=str, help='Input TSV file column to read')
    delete_entries_column_parser.add_argument('values', type=str, help='Values to delete')               
               
    args = parser.parse_args()
    
    if args.cmd == 'genomes':
        cmd_genomes(args.input, args.output, args.log)
    elif args.cmd == 'entrez_download':
        cmd_entrez_download(args.input, args.column, args.db, args.format, args.ext, args.out_dir)
    elif args.cmd == 'export_column':
        df_export_column(read_tsv_file(args.input), args.column, args.out_file)
    elif args.cmd == 'export_SRA':
        cmd_export_SRA(read_tsv_file(args.input), args.platform, args.out_file)
    elif args.cmd == 'check_duplicates':
        print(','.join(df_check_duplicates(read_tsv_file(args.input), args.column)))
    elif args.cmd == 'delete_entries':
        cmd_delete_entries(args.input, args.column, args.values.split(','))

        
if __name__ == "__main__":
    main()
