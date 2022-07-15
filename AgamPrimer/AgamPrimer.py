import pandas as pd 
import allel 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import patches
import malariagen_data

ag3 = malariagen_data.Ag3()


def prepare_gDNA_sequence(target_loc, amplicon_size_range, genome_seq, assay_name):
    """
    Extracts sequence of interest from genome sequence
    """
    # Set up range for the input sequence, we'll take the middle range of the amplicon size and add that either
    # side of the target SNP
    start = target_loc - np.mean(amplicon_size_range)
    end = target_loc + np.mean(amplicon_size_range)  
    # join array into be one character string, and store the positions of these sequences for later
    target_sequence = ''.join(genome_seq[start:end].compute().astype(str))
    gdna_pos = np.arange(start, end).astype(int) + 1
    print(f"The target sequence is {len(target_sequence)} bases long")

    # We need the target snp indices within the region of interest
    target_loc_primer3 = int(np.where(gdna_pos == target_loc)[0])
    target_loc_primer3 = [target_loc_primer3, 10]
    print(f"the target snp is {target_loc_primer3[0]} bp into our target sequence")
    
    seq_parameters = {
            'SEQUENCE_ID': assay_name,
            'SEQUENCE_TEMPLATE': target_sequence,
        'SEQUENCE_TARGET':target_loc_primer3,
        'GENOMIC_SEQUENCE_TARGET': target_loc
        }
    
    return(target_sequence, gdna_pos, seq_parameters)

def prepare_cDNA_sequence(transcript, gff, genome_seq, assay_name):
    """
    Extract exonic sequence for our transcript and record exon-exon junctions
    """
    #subset gff to your gene
    gff = gff.query(f"type == 'exon' & Parent == @transcript") 
    # Get fasta sequence for each of our exons, and remember gDNA position
    seq = dict()
    gdna_pos = dict()
    for i, exon in enumerate(zip(gff.start, gff.end)):
        seq[i] = ''.join(np.array(genome_seq)[exon[0]-1:exon[1]].astype(str))
        gdna_pos[i] = np.arange(exon[0]-1, exon[1])
    #concatenate exon FASTAs into one big sequence
    gdna_pos = np.concatenate(list(gdna_pos.values()))
    target_mRNA_seq = ''.join(seq.values())

    # Get list of exon junction positions
    exon_junctions = np.array(np.cumsum(gff.end - gff.start))[:-1]
    exon_sizes = np.array(gff.end - gff.start)[:-1]
    exon_junctions_pos = [ex + gff.iloc[i, 3] for i, ex in enumerate(exon_sizes)]
    print(f"Exon junctions for {transcript}:" , exon_junctions, exon_junctions_pos, "\n")
    seq_parameters = {
            'SEQUENCE_ID': assay_name,
            'SEQUENCE_TEMPLATE': target_mRNA_seq,
        'SEQUENCE_OVERLAP_JUNCTION_LIST':list(map(int, exon_junctions)),
        'TRANSCRIPT':transcript}

    return(target_mRNA_seq, list(map(int, exon_junctions)), gdna_pos, seq_parameters)    

def rev_complement(seq):
    BASES ='NRWSMBDACGTHVKSWY'
    return ''.join([BASES[-j] for j in [BASES.find(i) for i in seq][::-1]])

def consecutive(data, stepsize=1):
    arr = np.split(data, np.where(np.diff(data) != stepsize)[0]+1)
    arr = [[a.min(), a.max()] for a in arr]
    return (arr)

def get_primer_alt_frequencies(primer_df, gdna_pos, pair, sample_set, assay, contig):
  """
  Find the genomic locations of pairs of primers, and runs span_to_freq
  to get allele frequencies at those locations
  """
  primer_loc_fwd = primer_df.loc['PRIMER_LEFT', str(pair)][0]
  primer_size_fwd = primer_df.loc['PRIMER_LEFT', str(pair)][1]
  primer_loc_rev = primer_df.loc['PRIMER_RIGHT', str(pair)][0]+1
  primer_size_rev = primer_df.loc['PRIMER_RIGHT', str(pair)][1]
  
  freq_arr, ref_arr, pos_arr = get_primer_arrays(contig=contig, gdna_pos=gdna_pos, sample_set=sample_set, assaytype=assay)

  freq_fwd = freq_arr[primer_loc_fwd:primer_loc_fwd+primer_size_fwd]
  freq_rev = np.flip(freq_arr[primer_loc_rev-primer_size_rev:primer_loc_rev])
  ref_fwd = ref_arr[primer_loc_fwd:primer_loc_fwd+primer_size_fwd]
  ref_rev = ref_arr[primer_loc_rev-primer_size_rev:primer_loc_rev]
  ref_rev = np.array(list(rev_complement(''.join(ref_rev))), dtype=str)
  pos_fwd = pos_arr[primer_loc_fwd:primer_loc_fwd+primer_size_fwd]
  pos_rev = np.flip(pos_arr[primer_loc_rev-primer_size_rev:primer_loc_rev])

  fwd_df = pd.DataFrame({'position': pos_fwd, 'base':ref_fwd, 'alt_frequency':freq_fwd}) # Make dataframe for plotting
  rev_df = pd.DataFrame({'position': pos_rev, 'base':ref_rev, 'alt_frequency':freq_rev}) # Make dataframe for plotting
  fwd_df['base_pos'] = fwd_df['base'] + "_" + fwd_df['position'].astype(str)
  assert fwd_df.shape[0] == primer_size_fwd, "Wrong size primers"
  rev_df['base_pos'] = rev_df['base'] + "_" + rev_df['position'].astype(str) 
  assert rev_df.shape[0] == primer_size_rev, "Wrong size primers"
  return(fwd_df, rev_df)

def get_primer_arrays(contig, gdna_pos, sample_set, assaytype):

  if assaytype == 'gDNA':
    span_str=f'{contig}:{gdna_pos.min()}-{gdna_pos.max()}'                        
    ref_arr = ag3.genome_sequence(region=span_str).compute().astype(str)     # get dna sequence for span
    geno = ag3.snp_genotypes(region=span_str, sample_sets=sample_set).compute() # get genotypes
    freqs = allel.GenotypeArray(geno).count_alleles().to_frequencies()                                # calculate allele frequencies
    freq_arr = freqs[:, 1:].sum(axis=1)
    pos_arr = gdna_pos
  elif assaytype == 'qPCR':
    freq_arr = np.array([])
    ref_arr = np.array([])
    pos_arr = np.array([])
    exon_spans = np.array(consecutive(gdna_pos)) +1
    for span in exon_spans:
      span_str=f'{contig}:{span[0]}-{span[1]}'                        
      primer_ref_seq = ag3.genome_sequence(region=span_str).compute().astype(str)     # get dna sequence for span
      geno = ag3.snp_genotypes(region=span_str, sample_sets=sample_set).compute() # get genotypes
      freqs = allel.GenotypeArray(geno).count_alleles().to_frequencies()                                # calculate allele frequencies
      freq_arr = np.append(freq_arr, freqs[:, 1:].sum(axis=1))
      ref_arr = np.append(ref_arr, primer_ref_seq)
      pos_arr = np.append(pos_arr, np.arange(span[0], span[1]+1).astype(int))

  return(freq_arr, ref_arr, pos_arr)

def plot_pair_text(primer_df, pair, ax, side, dict_dfs):
  """
  Plot text relating to primer characteristics (Tm etc)
  """
  tm = np.round(primer_df.loc[f'PRIMER_{side}_TM', str(pair)], 2)       
  gc = np.round(primer_df.loc[f'PRIMER_{side}_GC_PERCENT', str(pair)], 2)
  span = f"{int(dict_dfs[pair]['position'].min())}-{int(dict_dfs[pair]['position'].max())}"
  # Write text to plot for Tm, GC, span, and 3/5'
  ax.text(x=2, y=0.5, s=f"TM={tm}", bbox=dict(facecolor='white', edgecolor='black', alpha=0.4, boxstyle='round', pad=0.5))
  ax.text(x=2, y=0.7, s=f"GC={gc}", bbox=dict(facecolor='white', edgecolor='black', alpha=0.4, boxstyle='round', pad=0.5))
  ax.text(x=2, y=0.9, s=span, bbox=dict(facecolor='white', edgecolor='black', alpha=0.4, boxstyle='round', pad=0.5))
  ax.text(x=0.5, y=0.9, s="5'")
  ax.text(x=18, y=0.9, s="3'")


def plot_primer(primer_df, ax, i, dict_dfs, assay, side='LEFT', exon_junctions=None):
  """
  Plot primer allele frequencies and text
  """
  sns.scatterplot(ax=ax, x=dict_dfs[i]['base_pos'], y=dict_dfs[i]['alt_frequency'], s=200) 
  ax.set_xticklabels(dict_dfs[i]['base'])
  ax.set_ylim(0,1)
  ax.set_xlabel("")
  ax.set_ylabel("Alternate allele frequency")
  if side == 'RIGHT': 
    ax.set_ylabel("")
    ax.set_yticklabels("")

  if assay == 'qPCR':
    a = primer_df.loc[f'PRIMER_{side}', str(i)]
    a = a[0] - a[1] if side == 'RIGHT' else a[0] + a[1]
    arr = np.array(exon_junctions)- a
    if (np.abs(arr) < 18).any():
      bases_in_new_exon  = arr[np.abs(arr) < 18]
      plt.setp(ax.get_xticklabels()[bases_in_new_exon[0]:], backgroundcolor="antiquewhite") if side == 'RIGHT' else plt.setp(ax.get_xticklabels()[:bases_in_new_exon[0]], backgroundcolor="antiquewhite")

  ax.set_title(f"FWD primer {i}") if side == 'LEFT' else  ax.set_title(f"REV primer {i}") 
  plot_pair_text(primer_df, i, ax, side, dict_dfs)


def plot_primer_ag3_frequencies(primer_df, gdna_pos, contig, sample_set, n_primer_pairs, assay_type, seq_parameters, save=True):
  """
  Loop through n primer pairs, retrieving frequency data and plot allele frequencies
  """
  name = seq_parameters['SEQUENCE_ID']
  exon_junctions = seq_parameters['SEQUENCE_OVERLAP_JUNCTION_LIST'] if assay_type == 'qPCR' else None
  transcript = seq_parameters['TRANSCRIPT'] if assay_type == 'qPCR' else None
  target_loc = seq_parameters['GENOMIC_SEQUENCE_TARGET'] if assay_type == 'gDNA' else None
  di_fwd = {}
  di_rev = {}
  # Loop through each primer pair and get the frequencies of alternate alleles, storing in dict
  for i in range(n_primer_pairs):
    di_fwd[i], di_rev[i] = get_primer_alt_frequencies(primer_df, gdna_pos, i, sample_set, assay_type, contig)

  # Plot data
  fig, ax = plt.subplots(n_primer_pairs , 2, figsize=[14, (n_primer_pairs*2)+2], constrained_layout=True)    
  if assay_type == 'gDNA':
    fig.suptitle(f"{name} primer pairs | {sample_set} | target {target_loc} bp", fontweight='bold')
  elif assay_type == 'qPCR':
    fig.suptitle(f"{name} primer pairs | {sample_set} | target {transcript}", fontweight='bold')

  for i in range(n_primer_pairs):
    plot_primer(primer_df, ax[i,0], i, di_fwd, assay_type, side='LEFT', exon_junctions=exon_junctions)
    plot_primer(primer_df, ax[i,1], i, di_rev, assay_type, side='RIGHT', exon_junctions=exon_junctions)
  if save: fig.savefig(f"{name}.{assay_type}.primers.png")
  return(di_fwd, di_rev)


def primer3_run_statistics(primer_dict):
  primer_df = pd.DataFrame.from_dict(primer_dict.items())          # Convert the dict into a pandas dataframe
  primer_df = primer_df.rename(columns={0:'parameter', 1:'value'}) # Rename the columns
  explanations_df = primer_df.iloc[:7, :]                          # Take the first 7 rows which are general
  for idx, row in explanations_df.iterrows():                      # Loop through each row and print information
      print(row['parameter'], " : ", row['value'], "\n")

def primer3_to_pandas(primer_dict):
  primer_df = pd.DataFrame.from_dict(primer_dict.items())          # Convert the dict into a pandas dataframe
  primer_df = primer_df.rename(columns={0:'parameter', 1:'value'}) # Rename the columns
  # Create a column which is primer pair #, and a column for primer parameter which does not contain primer pair #
  primer_df = primer_df.iloc[7:, :].copy()
  primer_df['primer_pair'] = primer_df['parameter'].str.extract("([0-9][0-9]|[0-9])")
  primer_df['parameter'] = primer_df['parameter'].str.replace("(_[0-9][0-9]|_[0-9])", "", regex=True)

  # Put the different primer pairs in different columns
  primer_df = primer_df.pivot(index='parameter', columns='primer_pair', values='value')

  # Get a list of the rows we need 
  primer_span = ['PRIMER_LEFT', 'PRIMER_RIGHT']
  required_info = ['SEQUENCE', 'TM', 'GC_PERCENT']
  required_info = [p + "_" + y for y in required_info for p in primer_span] + primer_span + ['PRIMER_PAIR_PRODUCT_SIZE']

  # Subset data frame
  primer_df = primer_df.loc[required_info, np.arange(primer_df.shape[1]).astype(str)]
  return(primer_df)

def get_gDNA_locs(gff, contig, start, end):
    locgff = gff.query("contig == @contig & type == 'exon' & start > @start & end < @end")
    min_= locgff.start.min() - 100
    max_ = locgff.end.max() + 100
    genegff = gff.query("contig == @contig & type == 'gene' & start > @min_ & end < @max_")
    return(locgff, min_, max_, genegff)

def get_qPCR_locs(gff, contig, transcript):
    # Load geneset (gff)
    locgff = gff.query("Parent == @transcript & type == 'exon'")
    min_= locgff.start.min() - 200
    max_ = locgff.end.max() + 200
    genegff = gff.query("contig == @contig & type == 'gene' & start > @min_ & end < @max_")
    return(locgff, min_, max_, genegff)

def plot_primer_locs(gff, contig, seq_parameters, n_primer_pairs, ax, di_fwd, di_rev, assay_type):
    # Load geneset (gff)
    if assay_type == 'gDNA':
      start = seq_parameters['GENOMIC_SEQUENCE_TARGET'] - 500
      end = seq_parameters['GENOMIC_SEQUENCE_TARGET'] + 500
      locgff, min_, max_, genegff = get_gDNA_locs(gff, contig, start, end)
    elif assay_type == 'qPCR':
      transcript = seq_parameters['TRANSCRIPT']
      locgff, min_, max_, genegff = get_qPCR_locs(gff, contig, transcript)
    # configure axes
    ax.set_xlim(min_, max_)
    ax.set_ylim(-0.5, 1.5)
    ax.ticklabel_format(useOffset=False)
    ax.axhline(0.5, color='k', linewidth=0.7, linestyle='--')
    sns.despine(ax=ax, left=True, bottom=False)
    #ax.set_yticks(ticks=[0.2,1.2], size=20)#labels=['- ', '+']
    ax.tick_params(top=False,left=False,right=False,labelleft=False,labelbottom=True)
    ax.tick_params(axis='x', which='major', labelsize=13)
    ax.set_ylabel("Exons")
    ax.set_xlabel(f"Chromosome {contig} position", fontdict={'fontsize':14})
    # Add rectangles for exons one at a time 
    for _, exon in locgff.iterrows():
        start, end = exon[['start', 'end']]
        e_name = exon['Name'][-2:]
        strand = exon['strand']
        if strand == '+':
            rect = patches.Rectangle((start, 0.55), end-start, 0.3, linewidth=3,
                                edgecolor='none', facecolor="grey", alpha=0.9)
            ax.text((start+end)/2, 0.65, e_name)
        else:
            rect = patches.Rectangle((start, 0.45), end-start, -0.3, linewidth=3,
                                edgecolor='none', facecolor="grey", alpha=0.9)
            ax.text((start+end)/2, 0.3, e_name)

        ax.add_patch(rect)
    for _, gene in genegff.iterrows():
        start, end = gene[['start', 'end']]
        size = end-start
        corr = size/4
        strand = gene['strand']
        if strand == '+':
            rect = patches.Rectangle((start, 0.55), end-start, 0.3, linewidth=3,
                                edgecolor='black', facecolor="none")
            ax.text(((start+end)/2)-corr, 0.95, s=gene['ID'], fontdict= {'fontsize':12}, weight='bold')
        else:
            rect = patches.Rectangle((start, 0.45), end-start, -0.3, linewidth=3,
                                edgecolor='black', facecolor="none")
            ax.text(((start+end)/2)-corr,  -0.3, s=gene['ID'], fontdict= {'fontsize':12},  weight='bold')
        ax.add_patch(rect)

    pal = sns.color_palette("Set2", n_primer_pairs)
    handles, labels = ax.get_legend_handles_labels()
    for pair in range(n_primer_pairs):
      lower_fwd, upper_fwd = di_fwd[pair]['position'].min() , di_fwd[pair]['position'].max()
      lower_rev, upper_rev = di_rev[pair]['position'].min() , di_rev[pair]['position'].max()

      plt.arrow(lower_fwd, 0.8+(2/(10-(pair))), upper_fwd-lower_fwd, 0, width=0.03, length_includes_head=True, color=pal[pair])
      plt.arrow(upper_rev, 0.8+(2/(10-(pair))), lower_rev-upper_rev, 0, width=0.03, length_includes_head=True, color=pal[pair])
      # manually define a new patch 
      patch = patches.Patch(color=pal[pair], label=f'pair {pair}')
      # handles is a list, so append manual patch
      handles.append(patch) 
    # plot the legend
    plt.legend(handles=handles, loc='best')
