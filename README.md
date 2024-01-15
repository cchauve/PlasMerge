# PlasMerge

PlasMerge, a tool to merge plasmid bins.

**TO DO** Better introduction, overview of problem and of the method

## Running PlasMerge

'''
python src/plasmerge.py \
       -s <sample name> \
       -a <gzipped GFA file> \
       -p <plasmid scores TSV file> \
       -b <plasmid bins TSV file> \
       -r <plasmid bin source> \
       -d <directory where Gurobi sol file is written> \
       -f <path to output TSV file> \
       -g <GC intervals file; optional> \
       -t <threshold; optional, default=0.05> \
       -o <plasmid score offset; optional, default=0.5>
'''

Parameters description:
- 'sample name': string (example: 'SAMN32247522')
- 'gzipped GFA file': gzipped assembly graph in GFA format (example: 'test/gfas/SAMN32247522.gfa.gz')
- 'plasmid scores TSV file': file recording for each contig its plasmid score (between 0 and 1)  
  Format: '<contig ID><TAB><plasmid score>'
  Example: 'test/scores/SAMN32247522.scores.tsv'
- 'plasmid bins TSV file': file recording the initial plasmid bins
  Format: '<plasmid id><TAB><comma separated list of contigs in format contig:muliplicity><TAB><copy_number>
  where 'copy_number' is optional and if absent the copy number of a plasmid bin is the minimum of the read depth
  of its contigs read in the GFA file.  
  Example: 'test/pls_bins/SAMN32247522.gp.tsv', 'test/pls_bins/SAMN32247522.pbf.tsv'
- 'plasmid bin source': string (in experiments: 'gt' for ground truth, 'gp' for gplas2, 'mob' for MOB-recon, 'pbf' for PlasBin-flow)
- 'GC intervals file': path to file describing the GC content ratio intervals  
  Format: one floating number per line, increasing from '0' to '1'
  Default: '0, 0.4,  0.45, 0.5, 0.55, 0.6, 1'
- 'threshold': **Question** what is this parameter?
- 'plasmid score offset': float in [0,1]: used in the objective function

