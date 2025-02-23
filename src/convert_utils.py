"""
Functions converting plasmid files into the standard plasmid bins format
Format:
- line 1: header "plasmid contigs <optional:copy_number>"
- plasmid line: <plasmid id>TAB<comma-spearated list of contig_id:multiplicity>TAB<copy number>
"""

import pandas as pd
from collections import defaultdict

def convert_pbf_ground_truth(in_pbf_file, out_pls_file):
    """
    Convert PlasBin-flow ground truth file in_pbf_file into standard file out_pls_file
    plasBin-flow ground truth format:
    plasmid contig  contig_len
    GT_2    21      78327
    GT_2    32      22479
    ...
    """
    plasmids_df = pd.read_csv(
        in_pbf_file, sep='\t', header=0, skip_blank_lines=True, index_col=False
    )
    result = defaultdict(lambda: defaultdict(int))
    for idx,row in plasmids_df.iterrows():
        pls_id = row['plasmid']
        ctg_id = row['contig']
        result[pls_id][ctg_id] += 1
    with open(out_pls_file, 'w') as out_file:
        out_file.write('plasmid\tcontigs\n')
        for pls_id,ctgs_dict in result.items():
            ctgs = ','.join(
                [
                    f'{ctg_id}:{ctg_mult}'
                    for ctg_id,ctg_mult in ctgs_dict.items()
                ]
            )
            out_file.write(f'{pls_id}\t{ctgs}\n')

def convert_pbf_output(in_pbf_file, out_pls_file):
    """
    Convert PlasBin-flow output file in_pbf_file into standard file out_pls_file
    plasBin-flow output format:
    #Pls_ID Flow    GC_bin          Contigs
    P1              0.89    0.5-0.55                42:1.0,183:3.0,180:2.0,161:3.0,...
    ...
    """
    plasmids_df = pd.read_csv(
        in_pbf_file, sep='\t', header=0, skip_blank_lines=True, index_col=False
    )
    with open(out_pls_file, 'w') as out_file:
        out_file.write('plasmid\tcontigs\tcopy_number\n')
        for idx,row in plasmids_df.iterrows():
            pls_id = row['#Pls_ID']
            copy_number = row['Flow']
            ctgs = row['Contigs']
            out_file.write(f'{pls_id}\t{ctgs}\t{copy_number}\n')

def convert_gplas_output(in_gplas_file, out_pls_file):
    """
    Convert gplas output file in_gplas_file into standard file out_pls_file
    gplas output format:
    number Contig_name Prob_Chromosome Prob_Plasmid Prediction length coverage Bin
    55 S55_LN:i:2123_dp:f:4.363318193712384 0.03 0.97 Repeat 2123 4.36 1
    ...
    """
    plasmids_df = pd.read_csv(
        in_gplas_file, sep=' ', header=0, skip_blank_lines=True
    )
    result = defaultdict(lambda: defaultdict(int))
    for idx,row in plasmids_df.iterrows():
        pls_id = row['Bin']
        ctg_id = row['Contig_name'].split('_LN')[0][1:]
        # Alternate to check: ctg_id = row['number']
        result[pls_id][ctg_id] += 1
    with open(out_pls_file, 'w') as out_file:
        out_file.write('plasmid\tcontigs\n')
        for pls_id,ctgs_dict in result.items():
            if pls_id != 'Unbinned':
                ctgs = ','.join(
                    [
                        f'{ctg_id}:{ctg_mult}'
                        for ctg_id,ctg_mult in ctgs_dict.items()
                    ]
                )
                out_file.write(f'{pls_id}\t{ctgs}\n')

def convert_mobsuite_output(in_mob_file, out_pls_file):
    """
    Convert MOB-suite output file in_mob_file into standard file out_pls_file
    MOB-suite output format:
    sample_id       molecule_type   primary_cluster_id      secondary_cluster_id    contig_id       size    gc      md5     circularity_status    rep_type(s)     rep_type_accession(s)   relaxase_type(s)        relaxase_type_accession(s)      mpf_type        mpf_type_accession(s)   orit_type(s)    orit_accession(s)       predicted_mobility      mash_nearest_neighbor   mash_neighbor_distance  mash_neighbor_identification    repetitive_dna_id       repetitive_dna_type     filtering_reason
    SAMN32247302.fna        chromosome      -       -       1       437536  0.49618317121334016     5a52cb988e024881417d16ef9e3b9a94        not tested      -       -       -       -       -       -       -       -       -       -       -       -       -       -       -
    SAMN32247302.fna        plasmid AA176   AH853   21      78327   0.528132061741162       69dc0707165ecceede88636dff1b1bad        not tested    IncFIC,rep_cluster_2244 AP001918,CP033091_00053 MOBF    NC_017627_00068 MPF_F,MPF_F,MPF_F,MPF_F,MPF_F,MPF_F,MPF_F,MPF_F,MPF_F,MPF_T,MPF_Unknown   NC_009837_00049,NC_010488_00021,NC_018966_00040,NC_017639_00100,NC_007675_00027,NC_017639_00094,NC_019094_00090,NC_010409_00124,NC_022651_00077,NC_013437_00116,08-5333_00200   MOBF    HG796403        -       MG825378        0.004991        Escherichia coli        -       -       -
    ...
    """
    plasmids_df = pd.read_csv(
        in_mob_file, sep='\t', header=0, skip_blank_lines=True
    )
    result = defaultdict(lambda: defaultdict(int))
    for idx,row in plasmids_df.iterrows():
        if row['molecule_type'] == 'plasmid':
            pls_id = row['primary_cluster_id']
            ctg_id = row['contig_id']
            result[pls_id][ctg_id] += 1
    with open(out_pls_file, 'w') as out_file:
        out_file.write('plasmid\tcontigs\n')
        for pls_id,ctgs_dict in result.items():
            ctgs = ','.join(
                [
                    f'{ctg_id}:{ctg_mult}'
                    for ctg_id,ctg_mult in ctgs_dict.items()
                ]
            )
            out_file.write(f'{pls_id}\t{ctgs}\n')


if __name__ == "__main__":
    import os
    import sys

    # Converting a PlasBi-flow file
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    convert_pbf_output(in_file, out_file)
    
    # samples = ['SAMN32247302', 'SAMN32247345', 'SAMN32247425', 'SAMN32247519', 'SAMN32247522']
    # root = os.path.normpath('../test')

    # for sample in samples:
    #     print(f'SAMPLE: {sample}')
        
    #     in_file = os.path.join(root, 'ground_truth', f'{sample}_gt.tsv')
    #     out_file = os.path.join(root, 'pls_bins', f'{sample}.gt.tsv')
    #     convert_pbf_ground_truth(in_file, out_file)

    #     in_file = os.path.join(root, 'pbf', f'{sample}.pred.txt')
    #     out_file = os.path.join(root, 'pls_bins', f'{sample}.pbf.tsv')
    #     convert_pbf_output(in_file, out_file)
        
    #     in_file = os.path.join(root, 'gplas', f'{sample}.gplas2.tab')
    #     out_file = os.path.join(root, 'pls_bins', f'{sample}.gp.tsv')
    #     convert_gplas_output(in_file, out_file)
        
    #     in_file = os.path.join(root, 'mob', f'{sample}.contig_report.txt')
    #     out_file = os.path.join(root, 'pls_bins', f'{sample}.mob.tsv')
    #     convert_mobsuite_output(in_file, out_file)
