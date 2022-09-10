import pandas as pd 
import allel 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import patches
import malariagen_data
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import gget
import primer3

ag3 = malariagen_data.Ag3(url='gs://vo_agam_release/')

def prepare_gDNA_sequence(target_loc, amplicon_size_range, genome_seq, assay_name, assay_type, probe_exclude_region_size=20):
    """
    Extracts sequence of interest from genome sequence
    """
    # Set up range for the input sequence, we'll take the max range of the amplicon size and add that either
    # side of the target SNP
    amp_size_range = int(np.max(amplicon_size_range)) 
    start = target_loc - amp_size_range
    end = target_loc + amp_size_range
    # join array into be one character string, and store the positions of these sequences for later
    target_sequence = ''.join(genome_seq[start:end-1].compute().astype(str))
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
        
    if 'probe' in assay_type:
        seq_parameters['SEQUENCE_INTERNAL_EXCLUDED_REGION'] = [[1,target_loc_primer3[0]-probe_exclude_region_size], [target_loc_primer3[0]+probe_exclude_region_size, len(target_sequence)-(target_loc_primer3[0]+probe_exclude_region_size)]]
    
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




def primer_params(assay_type, primer_parameters=None, n_primer_pairs=None, amplicon_size_range=None, generate_defaults=False):

  if generate_defaults:
    primer_parameters  =  {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_TASK':'generic',
        'PRIMER_MIN_SIZE': 17,
        'PRIMER_MAX_SIZE': 24,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 30.0,
        'PRIMER_MAX_GC': 70.0,
        'PRIMER_MIN_THREE_PRIME_DISTANCE':3,          # this parameter is the minimum distance between successive pairs. If 1, it means successive primer pairs could be identical bar one base shift
        'PRIMER_INTERNAL_OPT_SIZE': 16,               # Probe size preferences if selected, otherwise ignored
        'PRIMER_INTERNAL_MIN_SIZE': 10,
        'PRIMER_INTERNAL_MAX_SIZE': 22,
        'PRIMER_INTERNAL_MIN_TM': 45,
        'PRIMER_INTERNAL_MAX_TM':65,                # Probe considerations are quite relaxed, assumed that LNAs will be used later to affect TM
        # Extra primer3 parameters can go here
        # In the same format as above                       
    }

  primer_parameters['PRIMER_NUM_RETURN'] = n_primer_pairs
  primer_parameters['PRIMER_PRODUCT_SIZE_RANGE'] = amplicon_size_range

  if assay_type == 'gDNA primers + probe':
      primer_parameters['PRIMER_PICK_INTERNAL_OLIGO'] = 1
      primer_parameters['PRIMER_PICK_RIGHT_PRIMER'] = 1
      primer_parameters['PRIMER_PICK_LEFT_PRIMER'] = 1
  elif assay_type == 'probe':
      primer_parameters['PRIMER_PICK_INTERNAL_OLIGO'] = 1
      primer_parameters['PRIMER_PICK_RIGHT_PRIMER'] = 0
      primer_parameters['PRIMER_PICK_LEFT_PRIMER'] = 0
  return(primer_parameters)

def primer3_run_statistics(primer_dict, assay_type):
    _, row_start = return_oligo_list(assay_type)
    primer_dict = convert_results_dict_naming(primer_dict)
    primer_df = pd.DataFrame.from_dict(primer_dict.items())          # Convert the dict into a pandas dataframe
    primer_df = primer_df.rename(columns={0:'parameter', 1:'value'}) # Rename the columns
    explanations_df = primer_df.iloc[:row_start, :]         # Take the first 7 rows which are general
    for idx, row in explanations_df.iterrows():                      # Loop through each row and print information
        print(row['parameter'], " : ", row['value'], "\n")


def return_oligo_list(assay_type):
    if assay_type == 'probe':
        oligos = ['probe']
        row_start = 5
    elif any(item == assay_type for item in ['gDNA primers', 'qPCR primers']):
        oligos = ['forward', 'reverse']
        row_start = 7
    elif assay_type == 'gDNA primers + probe':
        oligos = ['forward', 'reverse', 'probe']
        row_start = 8
    return(oligos, row_start)

def convert_results_dict_naming(primer_dict):

    k = {}
    for key in primer_dict.keys():
        if 'LEFT' in key:
            nkey = key.replace("LEFT", 'forward')
        elif 'RIGHT' in key:
            nkey = key.replace("RIGHT", 'reverse')
        elif 'INTERNAL' in key:
            nkey = key.replace("INTERNAL", 'probe')
        else:
            nkey = key
        k[nkey.lower()] = primer_dict[key]
    return(k)



def primer3_to_pandas(primer_dict, assay_type):
    oligos, row_start = return_oligo_list(assay_type)

    primer_dict = convert_results_dict_naming(primer_dict)
    primer_df = pd.DataFrame.from_dict(primer_dict.items())          # Convert the dict into a pandas dataframe
    primer_df = primer_df.rename(columns={0:'parameter', 1:'value'}) # Rename the columns
    # Create a column which is primer pair #, and a column for primer parameter which does not contain primer pair #
    primer_df = primer_df.iloc[row_start:, :].copy()
    primer_df['primer_pair'] = primer_df['parameter'].str.extract("([0-9][0-9]|[0-9])")
    primer_df['parameter'] = primer_df['parameter'].str.replace("(_[0-9][0-9]|_[0-9])", "", regex=True)
    
    #start_row = primer_df.eval('parameter == "primer_pair_penalty"')
    #if any(start_row):
    #  row_start = [i for i, x in enumerate(start_row) if x][0]
    #  primer_df = primer_df.iloc[row_start:, :]
    # P#ut the different primer pairs in different columns
    primer_df = primer_df.pivot(index='parameter', columns='primer_pair', values='value')

    # Get a list of the rows we need 
    primer_span = [f'primer_{oligo}' for oligo in oligos]
    required_info = ['sequence', 'TM', 'GC_PERCENT']
    required_info = [p + "_" + y for y in required_info for p in primer_span] + primer_span
    required_info = required_info + ['primer_PAIR_PRODUCT_SIZE'] if assay_type != 'probe' else required_info
    required_info = [string.lower() for string in required_info]

    # Subset data frame
    primer_df = primer_df.loc[required_info, np.arange(primer_df.shape[1]).astype(str)]
    return(primer_df)










##### PLOTTING

def complement(x):
  if x == 'A':
    return("T")
  elif x == 'C':
    return("G")
  elif x == 'G':
    return("C")
  elif x == 'T':
    return("A")

complement = np.vectorize(complement)

def rev_complement(seq):
    BASES ='NRWSMBDACGTHVKSWY'
    return ''.join([BASES[-j] for j in [BASES.find(i) for i in seq][::-1]])

def consecutive(data, stepsize=1):
    arr = np.split(data, np.where(np.diff(data) != stepsize)[0]+1)
    arr = [[a.min(), a.max()] for a in arr]
    return (arr)
  
def addZeroCols(freqs):
  freqlength = freqs.shape[1]
  needed = 4-freqlength
  if needed > 0:
    for i in range(needed):
      freqs = np.column_stack([freqs, np.repeat(0, freqs.shape[0])])
  return(freqs)

def get_base_freqs(freqs, ref_alt_array):
  assert freqs.shape == ref_alt_array.shape, "Shape of arrays is different"
  freq_df = pd.DataFrame(ref_alt_array)
  for i_base in range(4):
    for i in range(freqs.shape[0]):
      base = ref_alt_array[i, i_base]
      freq_df.loc[i, f'{base}_freq'] = freqs[i, i_base]
  return(freq_df)

def get_primer_arrays(contig, gdna_pos, sample_set, assay_type, sample_query=None):

  if any(item in assay_type for item in ['gDNA', 'probe']):
    span_str=f'{contig}:{gdna_pos.min()}-{gdna_pos.max()}'                        
    snps = ag3.snp_calls(region=span_str, sample_sets=sample_set)     # get genotypes
    ref_alt_arr = snps['variant_allele'].compute().values
    geno = snps['call_genotype']
    freq_arr = allel.GenotypeArray(geno).count_alleles().to_frequencies() 
    pos_arr = gdna_pos
  elif assay_type == 'qPCR primers':
    freq_arr = []
    ref_alt_arr = []
    pos_arr = np.array([])
    exon_spans = np.array(consecutive(gdna_pos)) + 1
    for span in exon_spans:
      span_str = f'{contig}:{span[0]}-{span[1]}'                        
      snps = ag3.snp_calls(region=span_str, sample_sets=sample_set)     # get genotypes
      ref_alts = snps['variant_allele']
      geno = snps['call_genotype']
      freqs = allel.GenotypeArray(geno).count_alleles().to_frequencies()                                # calculate allele frequencies
      freqs = addZeroCols(freqs)
      freq_arr.append(freqs)
      ref_alt_arr.append(ref_alts)
      pos_arr = np.append(pos_arr, np.arange(span[0], span[1]+1).astype(int))
    freq_arr = np.concatenate(freq_arr)
    ref_alt_arr = np.concatenate(ref_alt_arr)

  return(freq_arr, ref_alt_arr.astype('U13'), pos_arr)



def get_primer_alt_frequencies(primer_df, gdna_pos, pair, sample_set, assay_type, contig, sample_query):
  """
  Find the genomic locations of pairs of primers, and runs span_to_freq
  to get allele frequencies at those locations
  """

  oligos, _ = return_oligo_list(assay_type)
  base_freqs, ref_alt_arr, pos_arr = get_primer_arrays(contig=contig, gdna_pos=gdna_pos, sample_set=sample_set, assay_type=assay_type, sample_query=sample_query)

  freq_arr = base_freqs[:, 1:].sum(axis=1)

  di = {}
  for oligo in oligos:
    primer_loc = primer_df.loc[f'primer_{oligo}', str(pair)][0]
    primer_loc = primer_loc + 1 if oligo == 'reverse' else primer_loc
    primer_size = primer_df.loc[f'primer_{oligo}', str(pair)][1]
    if oligo in ['forward', 'probe']:
      freq = freq_arr[primer_loc:primer_loc+primer_size]
      base_freqs_arr = base_freqs[primer_loc:primer_loc+primer_size, :]
      ref = ref_alt_arr[primer_loc:primer_loc+primer_size, 0]
      ref_alt = ref_alt_arr[primer_loc:primer_loc+primer_size, :]
      pos = pos_arr[primer_loc:primer_loc+primer_size]
    elif oligo == 'reverse':
      freq = np.flip(freq_arr[primer_loc-primer_size:primer_loc])
      base_freqs_arr = base_freqs[primer_loc-primer_size:primer_loc, :]
      base_freqs_arr = np.flip(base_freqs_arr, axis=0)
      ref = ref_alt_arr[primer_loc-primer_size:primer_loc, 0]
      ref = np.array(list(rev_complement(''.join(ref))), dtype=str)
      ref_alt = ref_alt_arr[primer_loc-primer_size:primer_loc, :]
      ref_alt = complement(np.flip(ref_alt, axis=0))
      pos = np.flip(pos_arr[primer_loc-primer_size:primer_loc])

    df = pd.DataFrame({'position': pos, 'base':ref, 'alt_frequency':freq}) # Make dataframe for plotting
    df['base_pos'] = df['base'] + "_" + df['position'].astype(str)
    assert df.shape[0] == primer_size, "Wrong size primers"

    freq_df = get_base_freqs(addZeroCols(base_freqs_arr), ref_alt).filter(like='freq')
    df = pd.concat([df, freq_df], axis=1)
    di[oligo] = df
  return(di)



def plotly_primers(primer_df, res_dict, name, assay_type, sample_set, target_loc, transcript, save=True):

  oligos, _ = return_oligo_list(assay_type)
  if len(oligos) == 2:
    plt_title = ['Forward primer', 'Reverse primer']
  elif len(oligos) == 3:
    plt_title = ['Forward primer', 'Reverse primer', 'Probe']
  elif len(oligos) == 1:
    plt_title = ['Probe']
    
  title_list = []
  for pair in primer_df:
      for oligo in plt_title:
        title_list.append(f"{oligo} {pair}")

  hover_template = "<br>".join(["Base / Position: %{customdata[4]}",
                                              "Total Alternate freq: %{y}",
                                              "A_freq: %{customdata[0]}",
                                              "C_freq: %{customdata[1]}",
                                              "G_freq: %{customdata[2]}",
                                              "T_freq: %{customdata[3]}"])
                              
  fig = make_subplots(rows=len(primer_df.columns), cols=len(oligos), subplot_titles=title_list, horizontal_spacing = 0.03, vertical_spacing=0.08)
  fig.update_annotations(font_size=13)
  for idx, oligo in enumerate(oligos):
    idx = idx+1
    for i in primer_df:
      i = int(i)
      row_i = i+1

      color = [-1 if v == 0 else 1 if v > 0 else 0 for v in res_dict[i][oligo]['alt_frequency']]
      colorscale = [[0, 'lightgray'], [0.5, 'lightgray'], [1, 'dodgerblue']]

      tm = np.round(primer_df.loc[f'primer_{oligo}_tm', str(i)], 2)       
      gc = np.round(primer_df.loc[f'primer_{oligo}_gc_percent', str(i)], 2)
      span = f"{int(res_dict[i][oligo]['position'].min())}-{int(res_dict[i][oligo]['position'].max())}"
      # Write text to plot for Tm, GC, span, and 3/5'

      fig.add_trace(go.Scatter(x=res_dict[i][oligo]['base_pos'], 
                              y=res_dict[i][oligo]['alt_frequency'], customdata=res_dict[i][oligo][['A_freq', 'C_freq', 'G_freq', 'T_freq', 'base_pos']], 
                              hovertemplate=hover_template,
                              mode='markers', marker=dict(size=14, color=color, colorscale=colorscale, line=dict(width=2, color='black')), marker_symbol='circle'), row=row_i, col=idx)
      fig.add_annotation(row=row_i, col=idx, x=res_dict[i][oligo]['base_pos'][0], y=0.8, text="5'", showarrow=False)
      fig.add_annotation(row=row_i, col=idx, x=res_dict[i][oligo]['base_pos'].to_numpy()[-1], y=0.8, text="3'", showarrow=False)
      fig.add_annotation(row=row_i, col=idx, x=res_dict[i][oligo]['base_pos'].to_numpy()[4], y=0.92, text=span, showarrow=False)
      fig.add_annotation(row=row_i, col=idx, x=res_dict[i][oligo]['base_pos'].to_numpy()[-7], y=0.92, text=f"GC={gc}", showarrow=False)
      fig.add_annotation(row=row_i, col=idx, x=res_dict[i][oligo]['base_pos'].to_numpy()[-3], y=0.92, text=f"TM={tm}", showarrow=False)
      
      fig.update_xaxes(row=row_i, col=idx, tickmode = 'array', tickvals = res_dict[i][oligo]['base_pos'], ticktext=res_dict[i][oligo]['base'], tickangle=0, mirror=True)
      if idx > 1:
        fig.update_yaxes(row=row_i, col=idx, range=[0,1], tickvals=np.arange(0, 1, 0.2), showticklabels=False, mirror=True)
      else:
        fig.update_yaxes(row=row_i, col=idx, tickvals=np.arange(0, 1, 0.2), range=[0,1], mirror=True)

      if ((row_i % 2) == 0) & (idx == 1):
        fig.update_yaxes(row=row_i, col=idx, title='Alternate allele frequency')

  if any(item in assay_type for item in ['gDNA']):
    title_text = f"{name} primer pairs | {sample_set} | target {target_loc} bp"
  elif assay_type == 'probe':
    title_text = f"{name} probe | {sample_set} | target {target_loc} bp"
  elif assay_type == 'qPCR primers':
      title_text = f"{name} primer pairs | {sample_set} | target {transcript}"

  #fig.update_traces(customdata=customdata, hovertemplate=hovertemplate)
  fig.update_layout(height=200*len(primer_df.columns), width=500*len(oligos),
                      title_text=title_text, title_x = 0.5, template="simple_white", showlegend=False)
  if save:
    fig.write_html(f"{name}_{assay_type}.html")
    fig.write_image(f"{name}_{assay_type}.pdf")
  fig.show()



def plot_primer_ag3_frequencies(primer_df, gdna_pos, contig, sample_set, assay_type, seq_parameters, save=True, sample_query=None):
  """
  Loop through n primer pairs, retrieving frequency data and plot allele frequencies
  """

  if sample_query != None:
    print(f"Subsetting allele frequencies to {sample_query}")

  name = seq_parameters['SEQUENCE_ID']
  exon_junctions = seq_parameters['SEQUENCE_OVERLAP_JUNCTION_LIST'] if assay_type == 'qPCR primers' else None
  transcript = seq_parameters['TRANSCRIPT'] if assay_type == 'qPCR primers' else None
  target_loc = seq_parameters['GENOMIC_SEQUENCE_TARGET'] if any(item in assay_type for item in ['gDNA', 'probe']) else None
  res_dict = {}
  # Loop through each primer pair and get the frequencies of alternate alleles, storing in dict
  for i in range(len(primer_df.columns)):
    res_dict[i] = get_primer_alt_frequencies(primer_df, gdna_pos, i, sample_set, assay_type, contig, sample_query)
  
  # Plot data with plotly
  plotly_primers(primer_df=primer_df, res_dict=res_dict, name=name, assay_type=assay_type, sample_set=sample_set, target_loc=target_loc, transcript=transcript, save=save)

  return(res_dict)








def get_gDNA_locs(gff, contig, start, end):
    locgff = gff.query("contig == @contig & type == 'exon' & start < @end & end > @start")
    min_= locgff.start.min() - 100
    max_ = locgff.end.max() + 100
    genegff = gff.query("contig == @contig & type == 'gene' & start < @end & end > @start")
    return(locgff, min_, max_, genegff)

def get_qPCR_locs(gff, contig, transcript):
    # Load geneset (gff)
    locgff = gff.query("Parent == @transcript & type == 'exon'")
    min_= locgff.start.min() - 200
    max_ = locgff.end.max() + 200
    genegff = gff.query("contig == @contig & type == 'gene' & start > @min_ & end < @max_")
    return(locgff, min_, max_, genegff)



def plot_primer_locs(primer_res_dict, primer_df, gff, contig, seq_parameters, assay_type, legend_loc='best', save=True):
    
    oligos, _ = return_oligo_list(assay_type)
    assay_name = seq_parameters['SEQUENCE_ID']
    # Load geneset (gff)
    if any(item in assay_type for item in ['gDNA', 'probe']):
        start = seq_parameters['GENOMIC_SEQUENCE_TARGET'] - 500
        end = seq_parameters['GENOMIC_SEQUENCE_TARGET'] + 500
        locgff, min_, max_, genegff = get_gDNA_locs(gff, contig, start, end)
        min_ = np.min([min_, start])
        max_ = np.max([max_, end])
    elif assay_type == 'qPCR primers':
        transcript = seq_parameters['TRANSCRIPT']
        locgff, min_, max_, genegff = get_qPCR_locs(gff, contig, transcript)

    if locgff.empty:
      print("No exons in close proximity for loc plot")
      return

    fig, ax = plt.subplots(1,1, figsize=[16,4])
    # configure axes
    if min_ in ['inf', 'NaN']:
      min_ = start
    if max_ in ['inf', 'NaN']:
      max_ = end

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
        ex_start, ex_end = exon[['start', 'end']]
        e_name = exon['Name'][-2:]
        strand = exon['strand']
        if strand == '+':
            rect = patches.Rectangle((ex_start, 0.55), ex_end-ex_start, 0.3, linewidth=3,
                                edgecolor='none', facecolor="grey", alpha=0.9)
            ax.text((ex_start+ex_end)/2, 0.65, e_name)
        else:
            rect = patches.Rectangle((ex_start, 0.45), ex_end-ex_start, -0.3, linewidth=3,
                                edgecolor='none', facecolor="grey", alpha=0.9)
            ax.text((ex_start+ex_end)/2, 0.3, e_name)
        ax.add_patch(rect)

    tot_genes = genegff.shape[0]
    for i, gene in genegff.reset_index(drop=True).iterrows():
        start, end = gene[['start', 'end']]
        diff = np.diff([min_, max_])
        interval = diff/tot_genes+1
        name_point = min_ + (interval*i+1)
        strand = gene['strand']
        if strand == '+':
            rect = patches.Rectangle((start, 0.55), end-start, 0.3, linewidth=3,
                                edgecolor='black', facecolor="none")
            ax.text(name_point, 0.95, s=gene['ID'], fontdict= {'fontsize':12}, weight='bold')
        else:
            rect = patches.Rectangle((start, 0.45), end-start, -0.3, linewidth=3,
                                edgecolor='black', facecolor="none")
            ax.text(name_point,  -0.3, s=gene['ID'], fontdict= {'fontsize':12},  weight='bold')
        ax.add_patch(rect)

    pal = sns.color_palette("Set2", len(primer_df.columns))
    handles, labels = ax.get_legend_handles_labels()
    for pair in primer_df:
      pair = int(pair)
      for oligo in oligos:
        lower, upper = primer_res_dict[pair][oligo]['position'].min() , primer_res_dict[pair][oligo]['position'].max()

        if oligo == 'forward':
            plt.arrow(lower, 0.8+(2/(10-(pair))), upper-lower, 0, width=0.03, length_includes_head=True, color=pal[pair])
        elif oligo == 'reverse':
            plt.arrow(upper, 0.8+(2/(10-(pair))), lower-upper, 0, width=0.03, length_includes_head=True, color=pal[pair])
        elif oligo == 'probe':
            ax.axhline(y=0.8+(2/(10-(pair))), xmin=lower, xmax=upper)
            line = plt.Line2D((lower, upper), (0.8+(2/(10-(pair))), 0.8+(2/(10-(pair)))), lw=2.5, color=pal[pair])
            ax.add_line(line)
        # manually define a new patch 
      patch = patches.Patch(color=pal[pair], label=f'pair {pair}')
      # handles is a list, so append manual patch
      handles.append(patch) 
    # plot the legend
    plt.legend(handles=handles, loc=legend_loc)
    if save:
      fig.savefig(f"{assay_name}_primer_locs.png", dpi=300)


# def batch_primer(target_df, assay_type, primer_parameters, amplicon_size_range, sample_sets, sample_query):
#     # Connect to the malariagen_data ag3 API
#     ag3 = malariagen_data.Ag3() #pre=True
    
#     genome_seq = {}
#     for idx, target in target_df.iterrows():
#         assay_name = target['name']
#         target_loc = target['target']
#         contig = target['contig']
#         genome_seq[contig] = ag3.genome_sequence(region=contig)
#         print(f"Our genome sequence for {contig} is {genome_seq[contig].shape[0]} bp long")
        
#         if any(item in assay_type for item in ['gDNA', 'probe']):
#           # genomic DNA
#           target_sequence, gdna_pos, seq_parameters = prepare_gDNA_sequence(target_loc=target_loc, amplicon_size_range=amplicon_size_range, genome_seq=genome_seq[contig], assay_name=assay_name, assay_type=assay_type)
#         elif assay_type == 'qPCR primers':
#           # RT-quantitative PCR, cDNA
#           target_sequence, exon_junctions, gdna_pos, seq_parameters = prepare_cDNA_sequence(transcript=target_loc, gff=ag3.geneset(), genome_seq=genome_seq, assay_name=assay_name)
        
#         primer_dict = primer3.designPrimers(seq_args=seq_parameters, global_args=primer_parameters)
#         primer_df = primer3_to_pandas(primer_dict, assay_type)
#         primer_df.to_csv(f"{assay_name}.{assay_type}.primers.tsv", sep="\t")
#         primer_df.to_excel(f"{assay_name}.{assay_type}.primers.xlsx")
        
#         results_dict = plot_primer_ag3_frequencies(primer_df=primer_df,
#                                                 gdna_pos=gdna_pos,
#                                                 contig=contig,
#                                                 sample_set=sample_sets, 
#                                                 sample_query=sample_query,
#                                                 assay_type=assay_type,
#                                                 seq_parameters=seq_parameters,
#                                                 save=True)
        
#         plot_primer_locs(primer_res_dict=results_dict, 
#                                     primer_df=primer_df, 
#                                     assay_type=assay_type, 
#                                     gff=ag3.geneset(), 
#                                     contig=contig, 
#                                     seq_parameters=seq_parameters, 
#                                     legend_loc='lower left',
#                                     save=True)





def gget_blat_genome(primer_df, assay_type, assembly='anoGam3'):
    oligos, _ = return_oligo_list(assay_type=assay_type)

    pair_dict = {}
    for primer_pair in primer_df:
        oligo_list = []
        for oligo in oligos:
            seq = primer_df[primer_pair].loc[f'primer_{oligo}_sequence']
            blat_df = gget.blat(sequence=seq, seqtype='DNA', assembly=assembly)
            if blat_df is None:
                print(f"No hit for {oligo} - pair {primer_pair}")
                continue
            blat_df.loc[:, 'primer'] = f"{oligo}_{primer_pair}"
            oligo_list.append(blat_df.set_index('primer'))
          
        if oligo_list:
          pair_dict[primer_pair] = pd.concat(oligo_list)
        elif not oligo_list:
          continue
    
    if pair_dict:
      return(pd.concat(pair_dict))
    else:
      print("No HITs found for these primer pairs")


def designPrimers(assay_type, assay_name, min_amplicon_size, max_amplicon_size, n_primer_pairs, contig, target, primer_parameters, sample_set, sample_query=None, save=True):

  amplicon_size_range = [[min_amplicon_size, max_amplicon_size]]
  ## adds some necessary parameters depending on assay type
  if primer_parameters == 'default':
    primer_parameters = primer_params(primer_parameters=None, assay_type=assay_type, n_primer_pairs=n_primer_pairs, amplicon_size_range=amplicon_size_range, generate_defaults=True)
  else:
    primer_parameters = primer_params(primer_parameters=primer_parameters, assay_type=assay_type, n_primer_pairs=n_primer_pairs, amplicon_size_range=amplicon_size_range, generate_defaults=False) 

  if assay_type == 'qPCR primers':
    assert not isinstance(target, int), "qPCR primers chosen but an AGAP identifier is not provided as the target"
    assert target.startswith("AGAP"), "qPCR primers chosen but an AGAP identifier is not provided as the target"
    transcript = target
    target_loc = ""
  else:
    assert isinstance(target, int), "For genomic DNA the target should be an integer within the contig"
    transcript = ""
    target_loc = target

  genome_seq = ag3.genome_sequence(region=contig)
  print(f"Our genome sequence for {contig} is {genome_seq.shape[0]} bp long")

  if any(item in assay_type for item in ['gDNA', 'probe']):
    # genomic DNA
    target_sequence, gdna_pos, seq_parameters = prepare_gDNA_sequence(target_loc=target_loc, amplicon_size_range=amplicon_size_range, genome_seq=genome_seq, assay_name=assay_name, assay_type=assay_type)
  elif assay_type == 'qPCR primers':
    # RT-quantitative PCR, cDNA
    target_sequence, exon_junctions, gdna_pos, seq_parameters = prepare_cDNA_sequence(transcript=transcript, gff=ag3.geneset(), genome_seq=genome_seq, assay_name=assay_name)
    
  primer_dict = primer3.designPrimers(seq_args=seq_parameters, global_args=primer_parameters)
  #AgamPrimer.primer3_run_statistics(primer_dict, assay_type)
  primer_df = primer3_to_pandas(primer_dict, assay_type)
  
  if save:
    primer_df.to_csv(f"{assay_name}.{assay_type}.primers.tsv", sep="\t")
    primer_df.to_excel(f"{assay_name}.{assay_type}.primers.xlsx")

  results_dict = plot_primer_ag3_frequencies(primer_df=primer_df,
                                              gdna_pos=gdna_pos,
                                              contig=contig,
                                              sample_set=sample_set, 
                                              sample_query=sample_query,
                                              assay_type=assay_type,
                                              seq_parameters=seq_parameters,
                                              save=save)
  plot_primer_locs(primer_res_dict=results_dict, primer_df=primer_df, assay_type=assay_type, gff=ag3.geneset(), contig=contig, seq_parameters=seq_parameters, legend_loc='lower left', save=save)
  blat_df = gget_blat_genome(primer_df, assay_type, assembly='anoGam3')
  return(primer_df, blat_df)











#def plot_primer(primer_df, ax, i, res_dict, assay_type, exon_junctions=None, target_loc=None):
#   """
#   Plot primer allele frequencies and text
#   """
#   oligos, _ = return_oligo_list(assay_type)

#   for idx, oligo in enumerate(oligos):
#     if len(oligos) > 1:
#       axes = ax[i, idx]
#     else:
#       axes = ax[i]

#     sns.scatterplot(ax=axes, x=res_dict[i][oligo]['base_pos'], y=res_dict[i][oligo]['alt_frequency'], s=200) 
#     axes.set_xticklabels(res_dict[i][oligo]['base'])
#     axes.set_ylim(0,1)
#     axes.set_xlabel("")
#     axes.set_ylabel("Alternate allele frequency")
#     if idx > 0: 
#       axes.set_ylabel("")
#       axes.set_yticklabels("")

#     if assay_type == 'qPCR primers':
#       a = primer_df.loc[f'primer_{oligo}', str(i)]
#       a = a[0] - a[1] if oligo == 'reverse' else a[0] + a[1]
#       arr = np.array(exon_junctions) - a
#       if (np.abs(arr) < 18).any():
#         bases_in_new_exon  = arr[np.abs(arr) < 18]
#         plt.setp(axes.get_xticklabels()[bases_in_new_exon[0]:], backgroundcolor="antiquewhite") if oligo == 'reverse' else plt.setp(axes.get_xticklabels()[:bases_in_new_exon[0]], backgroundcolor="antiquewhite")
    
#     if oligo == 'LEFT':
#       axes.set_title(f"forward primer {i}") 
#     elif oligo == 'reverse':
#       axes.set_title(f"reverse primer {i}")
#     elif oligo == 'INTERNAL':
#       axes.set_title(f"probe {i}")
#       idx = np.where(res_dict[i][oligo]['position'] == target_loc)[0][0]
#       plt.setp(axes.get_xticklabels()[idx], backgroundcolor="antiquewhite")
#     plot_pair_text(primer_df, i, axes, oligo, res_dict)


# def plot_pair_text(primer_df, pair, ax, oligo, res_dict):
#   """
#   Plot text relating to primer characteristics (Tm etc)
#   """
#   tm = np.round(primer_df.loc[f'primer_{oligo}_TM', str(pair)], 2)       
#   gc = np.round(primer_df.loc[f'primer_{oligo}_GC_PERCENT', str(pair)], 2)
#   span = f"{int(res_dict[pair][oligo]['position'].min())}-{int(res_dict[pair][oligo]['position'].max())}"
#   # Write text to plot for Tm, GC, span, and 3/5'
#   ax.text(x=2, y=0.5, s=f"TM={tm}", bbox=dict(facecolor='white', edgecolor='black', alpha=0.4, boxstyle='round', pad=0.5))
#   ax.text(x=2, y=0.7, s=f"GC={gc}", bbox=dict(facecolor='white', edgecolor='black', alpha=0.4, boxstyle='round', pad=0.5))
#   ax.text(x=2, y=0.9, s=span, bbox=dict(facecolor='white', edgecolor='black', alpha=0.4, boxstyle='round', pad=0.5))
#   ax.text(x=0.5, y=0.9, s="5'")
#   ax.text(x=18, y=0.9, s="3'")