#!./venv/bin/python

import os
import sys
import re
import subprocess
from io import StringIO
import weblogo #./venv/bin/pip install weblogo
import logomaker #./venv/bin/pip install logomaker
from Bio import pairwise2, SeqIO #./venv/bin/pip install biopython
from Bio.Align import substitution_matrices
# from Bio.Blast import NCBIXML
from Bio.Seq import reverse_complement
import numpy as np
import matplotlib.pyplot as plt
import requests
import pandas as pd
from scipy.stats import hypergeom, fisher_exact
import seaborn as sns
import copy


sample_datasets = {}

sample_datasets['PSF370_A1_v0'] = {}
sample_datasets['PSF370_A1_v0']['moise_folder'] = 'n/a'
sample_datasets['PSF370_A1_v0']['emote_barcode_file'] = './emote_barcodes/EmoteBarcode.1.csv'
sample_datasets['PSF370_A1_v0']['length_of_the_read'] = 50
sample_datasets['PSF370_A1_v0']['recognition_seq'] = 'CGGCACCAACCGAGG'
sample_datasets['PSF370_A1_v0']['recognition_mm'] = 0
sample_datasets['PSF370_A1_v0']['quantification_length'] = 7
sample_datasets['PSF370_A1_v0']['bases_used_in_quantification_seq'] = ['A', 'C', 'G']
sample_datasets['PSF370_A1_v0']['control_sequence'] = 'ACAGA'
sample_datasets['PSF370_A1_v0']['length_of_mapping_sequence'] = 20
sample_datasets['PSF370_A1_v0']['name_of_the_fastq_file'] = '../../raw/mpg_L31754-1_nEMOTE-Phage-qw_S1_R1_001.fastq.gz'
sample_datasets['PSF370_A1_v0']['name_of_the_output_file'] = '../../mpg_L31799-1_nEMOTE-Phage-R0_S1_R1_001'
sample_datasets['PSF370_A1_v0']['name_of_the_input_file'] = '../../mpg_L31799-1_nEMOTE-Phage-R0_S1_R1_001'
sample_datasets['PSF370_A1_v0']['name_of_the_output_files'] = '../../New_data_SF370_A1.fastq'
sample_datasets['PSF370_A1_v0']['name_of_genome_file'] = 'references/SF370_phage_a1_hybrid.fas'
sample_datasets['PSF370_A1_v0']['name_of_gene_information_file'] = 'references/SF370_phage_a1_hybrid.gff'
sample_datasets['PSF370_A1_v0']['name_of_sam_file'] = '../../Output.sam'
sample_datasets['PSF370_A1_v0']['kegg_org'] = 'spy'
sample_datasets['PSF370_A1_v0']['genbank_file'] = 'references/SF370_phage_a1_hybrid.gb'

class sample_data:
    def __init__(self):
        self.moise_folder = 'n/a'
        self.emote_barcode_file = './emote_barcodes/EmoteBarcode.1.csv'
        
        #InputParamsStepI:
        self.length_of_the_read = 50
        self.recognition_seq = 'CGGCACCAACCGAGG'
        self.recognition_mm = 0
        self.quantification_length = 7
        self.bases_used_in_quantification_seq = ['A', 'C', 'G']
        self.control_sequence = 'ACAGA'
        self.length_of_mapping_sequence = 20
        self.name_of_the_fastq_file = 'n/a'
        self.name_of_the_output_file = 'n/a'
        
        #InputParamsStepII:
        self.name_of_the_input_file = 'n/a'
        self.name_of_the_output_files = 'n/a'
        self.name_of_genome_file = 'n/a'
        self.name_of_gene_information_file = 'n/a'
        self.name_of_sam_file = 'n/a'
        
        self.kegg_org = 'n/a'
        self.genbank_file = 'n/a'

class Defaults:
      
    log_file = 'log.txt'
    
    cores = 10
    recognition_mm = 0
    pathway_enrich_adjust_p_threshold = 0.5
    motif = 'AAAAT'#'ARWAT'#cut site motif
    min_motif_per_CDS = 3
    
    min_end_value_to_report = 3#values below 3 will be set to 0
    
    # motif2 = '.GGC|C.GC|CG.C|CGG.'
    
    class programs:
        bowtie_build = '/array-1/tools/bowtie/1.1.2/bin/bowtie-build'
        bowtie_map = '/array-1/tools/bowtie/1.1.2/bin/bowtie'
        emote1_cmd = '/usr/bin/perl -Iscripts/ -MemoteInput -MemoteStat -McheckEmotes -MdetTimeMem ./scripts/emoteStepI.pl'
        emote2_py_cmd = './scripts/emote_step_II.py'
    
    # Definitions = {"D6A" : ("Empty plasmid", "Rep1", "-PNK"), "D6B" : ("Empty plasmid", "Rep1", "+PNK"), "D6C" : ("Empty plasmid", "Rep2", "-PNK"), "D6D" : ("Empty plasmid", "Rep2", "+PNK"), "D6E" : ("Empty plasmid", "Rep3", "-PNK"), "D6F" : ("Empty plasmid", "Rep3", "+PNK"), "D6H" : ("Toxin", "Rep1", "-PNK"), "D6I" : ("Toxin", "Rep1", "+PNK"), "D6J" : ("Toxin", "Rep2", "-PNK"), "D6K" : ("Toxin", "Rep2", "+PNK"), "D6L" : ("Toxin", "Rep3", "-PNK"), "D6M" : ("Toxin", "Rep3", "+PNK")}
    Definitions = {}
    Definitions_file = './references/replicate_definition.tsv'
    
    class logo_params:
        topx_list = [1000, 500]
        region_sizes = [15, 30]
        plot_types = ['png',]# 'eps']
        class set1:
            stacks_per_line = 61
            resolution = 300
            stack_aspect_ratio = 10
            yaxis_scale = 2
        class set2:
            stacks_per_line = 100
            resolution = 300
            stack_aspect_ratio = 10
            yaxis_scale = 1
    
    class Pairs:
        pass

def read_replicate_definition(datei, sampleID):
    Definition = {}
    header = []
    temp = {}
    sampleID = sampleID[1:]
    with open(datei, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            if header == []: header = line
            elif line[0].lower() == sampleID.lower():
                for e in range(len(line)):
                    if e>0:
                        if not header[e] in temp: temp[header[e]] = []
                        temp[header[e]].append(line[e])
    translate = {'empty_plasmid_noPNK': ("Empty plasmid", "-PNK"), 'empty_plasmid_PNK': ("Empty plasmid", "+PNK"), 'toxin_noPNK': ("Toxin", "-PNK"), 'toxin_PNK': ("Toxin", "+PNK")}
    
    for groupID in temp:
        for e,repID in enumerate(sorted(temp[groupID])):
            Definition[repID] = (translate[groupID][0], f"Rep{e+1}", translate[groupID][1])

    return Definition

def ambiguity_to_regex(string):
    translate = {'R':'[AG]', 'Y':'[CT]', 'S':'[GC]', 'W':'[AT]', 'K':'[GT]', 'M':'[AC]', 'B':'[CGT]', 'D':'[AGT]', 'H':'[ACT]', 'V':'[ACG]', 'N':'[ACGT]'}
    re_exp = []
    for letter in string:
        if letter in translate:
            re_exp.append(translate[letter])
        else:
            re_exp.append(letter)
    return re.compile(''.join(re_exp))

def plot_logo(data, to_type='information', ignore='', title='', ylim=[0, 1], background=[0.305, 0.195, 0.195, 0.305], colors='classic', vpad=0, figsize=[25, 2], filename='logo.png'):
    '''Code from Michael'''
    data_logo = logomaker.alignment_to_matrix(data, to_type=to_type, characters_to_ignore=ignore, background=background)
    data_logo.index += 1
    title = f"{title} (n = {len(data)})"
    logo = logomaker.Logo(data_logo, color_scheme=colors, vpad=vpad, figsize=figsize, fade_below=0, flip_below=False)
    max_bits = logo.df.sum(axis=1).max()
    logo.ax.set_ylim(0, max_bits*1.25)
    # logo.ax.set_ylim(ylim)
    logo.ax.set_title(title)
    logo.ax.set_xticks(np.arange(2, len(data[0])+1, 2))
    for typ in Defaults.logo_params.plot_types:
        temp_filename = filename
        if not temp_filename.endswith(f'.{typ}'): temp_filename += f'.{typ}'
        plt.savefig(temp_filename, dpi=300, format=typ, bbox_inches='tight')
    return temp_filename

def reverse_complement(string):
  table = str.maketrans('ACGTNWSMKRYBDHV.acgtnwsmkrybdhv','TGCANWSKMYRVHDBNtgcanwskmyrvhdb')
  string = string.translate(table)[::-1]
  return string

def create_reference_index(target=None):
    if target == None:
        for sampleID in dir(Defaults.Pairs):
            if not sampleID.startswith('__'):
                ref = getattr(Defaults.Pairs, sampleID).name_of_genome_file
                if not os.path.isfile(ref+'.1.ebwt'):
                    print(f'Creating index for {ref}')
                    p = subprocess.run([Defaults.programs.bowtie_build, ref, ref])
                    print('Error:', p.stderr)
    else:
        if not os.path.isfile(target+'.1.ebwt'):
            print(f'Creating index for {target}')
            p = subprocess.run([Defaults.programs.bowtie_build, target, target])
            print('Error:', p.stderr)
    return True

def run_emote_step1(sampleID):

    print(F'Running emote Step 1 for "{sampleID[1:]}"')

    if not os.path.exists('results'):
        os.makedirs('results')
    
    #create output folder
    output_folder = os.path.join('./results', sampleID[1:])
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    #create uncompressed copy of input fastq file (Emote cannot handle compressed input) and correcting the path in Defaults
    if not os.path.exists(os.path.join(output_folder, getattr(Defaults.Pairs, sampleID).name_of_the_fastq_file.split('/')[-1].replace('.gz',''))):
        cmd = f"zcat {getattr(Defaults.Pairs, sampleID).name_of_the_fastq_file} > {os.path.join(output_folder, getattr(Defaults.Pairs, sampleID).name_of_the_fastq_file.split('/')[-1].replace('.gz',''))}"
        print(cmd)
        p = subprocess.run(cmd, shell=True)
    getattr(Defaults.Pairs, sampleID).name_of_the_fastq_file = os.path.join(output_folder, getattr(Defaults.Pairs, sampleID).name_of_the_fastq_file.split('/')[-1].replace('.gz',''))
    
    #correcting output folder for step I in Defaults and other values
    getattr(Defaults.Pairs, sampleID).name_of_the_output_file = os.path.join(output_folder, 'emote_step1')
    getattr(Defaults.Pairs, sampleID).recognition_mm = Defaults.recognition_mm
    
    #create input csv file
    input_csv_file = os.path.join(output_folder, 'InputParamsStepI.csv')
    with open(input_csv_file, 'w') as outfile:
        lines = [[],[]]
        for param in ['length_of_the_read','recognition_seq','recognition_mm','quantification_length','bases_used_in_quantification_seq','control_sequence','length_of_mapping_sequence','name_of_the_fastq_file','name_of_the_output_file']:
            lines[0].append(param.replace('_', ' '))
            val = getattr(getattr(Defaults.Pairs, sampleID), param)
            if type(val) is list:
                val = ','.join(val)
            lines[1].append(str(val))
        for line in lines:
            outfile.write('\t'.join(line)+'\n')
    
    #run emote step 1
    if not os.path.exists(getattr(Defaults.Pairs, sampleID).name_of_the_output_file+'Emote.fastq'):
        cmd = f"{Defaults.programs.emote1_cmd} {getattr(Defaults.Pairs, sampleID).emote_barcode_file} {input_csv_file}"
        p = subprocess.run(cmd, shell=True)
        if p.stderr != None: print('Error:', p.stderr)
    
    return True

def run_bowtie_mapping(sampleID, mapping_referenceID=None):
    print(F'Running bowtie mapping for "{sampleID[1:]}"', flush=True)
    getattr(Defaults.Pairs, sampleID).name_of_sam_file = getattr(Defaults.Pairs, sampleID).name_of_the_output_file+'Emote.sam'
    if mapping_referenceID == None:
        mapping_referenceID = sampleID
    cmd = f"{Defaults.programs.bowtie_map} {getattr(Defaults.Pairs, mapping_referenceID).name_of_genome_file} {getattr(Defaults.Pairs, sampleID).name_of_the_output_file}Emote.fastq -S {getattr(Defaults.Pairs, sampleID).name_of_sam_file} -a --best --strata -v 1 -p {Defaults.cores}"
    print(cmd)
    if not os.path.exists(getattr(Defaults.Pairs, sampleID).name_of_sam_file):
        p = subprocess.run(cmd, shell=True)
        if p.stderr != None: print('Error:', p.stderr)
    
    dataset = ''
    with open(getattr(Defaults.Pairs, sampleID).name_of_sam_file, 'r') as infile:
        for line in infile:
            dataset += line.strip()
            if len(dataset) > 1000:
                break
    if len(dataset) == 0:
        print('Error: No data in sam file')
        print(getattr(Defaults.Pairs, sampleID).name_of_sam_file, '\nCmd:')
        print(cmd)
        sys.exit(1)
        return False
    
    return True

def run_emote_step2(sampleID):
    print(F'Running emote Step 2 for "{sampleID[1:]}"', flush=True)
    output_folder = os.path.join('./results', sampleID[1:])
    
    #adjusting parameter values
    getattr(Defaults.Pairs, sampleID).name_of_the_input_file = getattr(Defaults.Pairs, sampleID).name_of_the_output_file
    getattr(Defaults.Pairs, sampleID).name_of_the_output_files = os.path.join(output_folder, f'{sampleID[1:]}.new_data.')
    
    # #create input csv file
    # input_csv_file = os.path.join(output_folder, 'InputParamsStepII.csv')
    # with open(input_csv_file, 'w') as outfile:
    #     lines = [[],[]]
    #     for param in ['name_of_the_input_file','name_of_the_output_files','name_of_genome_file','name_of_gene_information_file','name_of_sam_file']:
    #         lines[0].append(param.replace('_', ' '))
    #         val = getattr(getattr(Defaults.Pairs, sampleID), param)
    #         if type(val) is list:
    #             val = ','.join(val)
    #         lines[1].append(str(val))
    #     for line in lines:
    #         outfile.write('\t'.join(line)+'\n')
    
    cmd = f"{Defaults.programs.emote2_py_cmd} --fq {getattr(Defaults.Pairs, sampleID).name_of_the_input_file}Emote.fastq -s {getattr(Defaults.Pairs, sampleID).name_of_sam_file} -r {getattr(Defaults.Pairs, sampleID).name_of_genome_file} -a {getattr(Defaults.Pairs, sampleID).name_of_gene_information_file} -b {getattr(Defaults.Pairs, sampleID).emote_barcode_file} -o {getattr(Defaults.Pairs, sampleID).name_of_the_output_files}   "
    # cmd = f"{Defaults.programs.emote2_cmd} {getattr(Defaults.Pairs, sampleID).emote_barcode_file} {input_csv_file}"
    # if not os.path.exists(getattr(getattr(Defaults.Pairs, sampleID), 'name_of_the_output_files')+'EmoteBarcodesReport.csv'):
    p = subprocess.run(cmd, shell=True)
    if p.stderr != None: print('Error:', p.stderr)
    return True

def remove_below_treshold_values_in_place(datei, threshold):
    dataset = []
    header = False
    cols = False
    count = 0
    total = 0
    with open(datei, 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            if not header:
                header = line
                cols = [e for e,c in enumerate(line) if c.startswith('D')]
                dataset.append('\t'.join(line))
            else:
                total += 1
                keep = False
                for p in cols:
                    if int(line[p]) < threshold:
                        line[p] = '0'
                    else:
                        keep = True
                if keep:
                    dataset.append('\t'.join(line))
                else:
                    count += 1
    
    with open(datei, 'w') as outfile:
        outfile.write('\n'.join(dataset)+'\n')

    # print(F"Removed {count:,} entries {(count*100/total):.2f}% due to min cov filter from {datei}")
    return datei

def run_follow_up_analysis(sampleID):
    print(F'Running follow-up analysis for "{sampleID[1:]}"', flush=True)
    
    #strand specific results
    pos_file = getattr(getattr(Defaults.Pairs, sampleID), 'name_of_the_output_files')+'UnambPosTable.csv'
    neg_file = getattr(getattr(Defaults.Pairs, sampleID), 'name_of_the_output_files')+'UnambNegTable.csv'

    remove_duplicate_entries_in_place(pos_file)
    remove_duplicate_entries_in_place(neg_file)
    
    remove_below_treshold_values_in_place(pos_file, Defaults.min_end_value_to_report)
    remove_below_treshold_values_in_place(neg_file, Defaults.min_end_value_to_report)

    #creating data frames and reading in data
    df = [pd.read_csv(pos_file, sep='\t'), pd.read_csv(neg_file, sep='\t')]#pos, neg
    
    
    #removing unwanted columns
    for i in range(2):
        df[i] = df[i][df[i].columns.difference(df[i].filter(regex='D6[A-Z]$').columns)]
    
    
    ####odd b/c this is missing in the analysis, but Moise tries to process it anyway and creates errors
    # df_totalQ_comb = pd.concat(df)
    # df_totalQ_comb = df_totalQ_comb[df_totalQ_comb.columns.difference(df_totalQ_comb.filter(regex='D6[A-Z]$').columns)]
    
    
    #combining replicates
    groups = {}
    for key, ID in [((l[0],l[2]),k) for k,l in Defaults.Definitions.items()]:
        #adding _Q to IDs to match input tables
        ID += '_Q'
        if not key in groups: groups[key] = []
        groups[key].append(ID)
        
    df_comb = []#pos, neg
    for i in range(2):
        df_comb.append(df[i][df[i].columns.difference(df[i].filter(regex='D6[A-Z]_Q').columns)].copy())
        for g in groups.keys():
            temp_group = [sub_g for sub_g in groups[g] if sub_g in df[i].columns]
            df_comb[i][g] = df[i][temp_group].mean(axis=1)

    df_totalQ_comb = pd.concat(df_comb, ignore_index=True)
    
    #filtering totals
    df_filtered = df_totalQ_comb.loc[(df_totalQ_comb[('Empty plasmid', '-PNK')] == 0) & (df_totalQ_comb[('Empty plasmid','+PNK')] == 0) & (df_totalQ_comb[('Toxin','-PNK')] == 0) & (df_totalQ_comb[('Toxin','+PNK')] >= 1)]
    
    #sorting and reordering columns
    df_sorted = df_filtered.sort_values(by=('Toxin','+PNK'), ascending=False)
    df_sorted = df_sorted.set_index('Genomic_pos')

    count_cols = df_sorted.filter(regex="PNK")
    df_sorted = df_sorted.drop(columns=count_cols.columns)
    for i,c in enumerate(count_cols.columns):
        df_sorted.insert(loc=i, column=c, value=count_cols[c])

    df_sorted.to_csv(f'results/{sampleID[1:]}/{sampleID[1:]}.df_sorted.csv',index=True, header=True)

    df_filtered = []#pos, neg
    df_sorted = []#pos, neg
    
    ## Filter-POSITIVE
    #in Moise's reference code, the first 3 columns were evaluated to '>=', apparently by mistake. They should evaluate to '== 0' as in the negative filter
    df_filtered.append(df_comb[0].loc[(df_comb[0][('Empty plasmid','-PNK')] == 0) & (df_comb[0][('Empty plasmid','+PNK')] == 0) & (df_comb[0][('Toxin','-PNK')] == 0) & (df_comb[0][('Toxin','+PNK')] >= 1)])

    ## Filter-NEGATIVE
    df_filtered.append(df_comb[1].loc[(df_comb[1][('Empty plasmid','-PNK')] == 0) & (df_comb[1][('Empty plasmid','+PNK')] == 0) & (df_comb[1][('Toxin','-PNK')] == 0) & (df_comb[1][('Toxin','+PNK')] >= 1)])

    for index in range(2):

        ## Sort by value of toxin +PNK
        df_sorted.append(df_filtered[index].sort_values(by=('Toxin','+PNK'), ascending=False).set_index('Genomic_pos'))

        ## Re-order columns to get the count columns to the front
        count_cols = df_sorted[index].filter(regex="PNK")
        df_sorted[index] = df_sorted[index].drop(columns=count_cols.columns)
        for i,c in enumerate(count_cols.columns):
            df_sorted[index].insert(loc=i, column=c, value=count_cols[c])

        suffix = 'pos' if index == 0 else 'neg'
        df_sorted[index].to_csv(f'results/{sampleID[1:]}/{sampleID[1:]}.target_{suffix}.csv',index=True, header=True)
    
    return df_sorted

def read_fasta(file):#returns [(ID, seq), ...]
    print(F"Reading fasta file: {file}")
    seqs = []
    with open(file, 'r') as infile:
        genome = ''
        title = ''
        for line in infile:
            if line.startswith('>') and len(genome) > 0:
                seqs.append((title, genome))
                genome = ''
                title = ''.join(line.strip()[1:].split())
            elif line.startswith('>'):
                title = ''.join(line.strip()[1:].split())
            else:
                genome += line.strip().upper()
        seqs.append((title, genome))
    return seqs

def remove_duplicate_entries_in_place(datei):
    dataset = set()
    header = False
    count_pre = 0
    with open(datei, 'r') as infile:
        for line in infile:
            line = tuple(line.strip().split('\t'))
            if not header:
                header = line
            else:
                count_pre += 1
                dataset.add((int(line[0]),line))
    dataset = [t[1] for t in sorted(dataset)]

    count_post = 0
    with open(datei, 'w') as outfile:
        outfile.write('\t'.join(header)+'\n')
        for line in dataset:
            count_post += 1
            outfile.write('\t'.join(line)+'\n')
    print(F"Removed {count_pre - count_post:,} duplicate entries from {datei}")
    return datei

def run_motif_scan(genome_seq, df_sorted, sampleID):
    for i in range(2):# pos, neg
        for topx in Defaults.logo_params.topx_list:
            for region_size in Defaults.logo_params.region_sizes:
                for s,logo_param_set in enumerate([Defaults.logo_params.set1, Defaults.logo_params.set2]):
                    seqs = []
                    try:
                        positions = df_sorted[i].iloc[0:topx, :].index
                        for pos in positions:
                            start = pos-region_size
                            end = pos+region_size
                            if start<0: start = 0
                            if end>len(genome_seq): end = len(genome_seq)
                            
                            for sub in df_sorted[i].loc[pos, "Strand"]:
                                if sub == 'minus':
                                    seqs.append(reverse_complement(genome_seq[start:end+1]).replace('T','U'))#end+1, b/c we want to have the cleavage position at the center of the plot
                                else:
                                    seqs.append(genome_seq[start:end+1].replace('T','U'))#end+1, b/c we want to have the cleavage position at the center of the plot
                                    break#b/c we only want to record this once per site
        
                        suffix = 'pos' if i == 0 else 'neg'
                        filename = F"results/{sampleID[1:]}/{sampleID[1:]}.top{topx}_region{region_size}_set{s}.{suffix}"
                        filename = plot_logo(seqs, filename=filename)
                        with open(filename+'.fas', 'w') as outfile:
                            for e, seq in enumerate(seqs):
                                outfile.write(F">{sampleID[1:]}_top{topx}_region{region_size}_set{s}_{suffix}_{e+1}\n{seq}\n")
                    except:
                        print(F"Unable to create logo for {sampleID[1:]} with top {topx}, region {region_size}, param set {s} for strand {suffix}; # of seqs: {len(seqs):,}", flush=True)
    return True

def run_codon_nucleotide_position(df_sorted, genome_seq, sampleID):
    genes_ds = []
    spy_gene_sequences = []
    for i in range(2):#pos, neg
        suffix = 'pos' if i == 0 else 'neg'
        positions_in_ORF = [] 
        positions_in_triplet =[]
        spy_gene_sequences.append({})
        genes = df_sorted[i].loc[pd.notna(df_sorted[i]['GeneSymbol|startPos|endPos|strand'])]
        genes_ds.append(genes)
        for pos,row in genes.iterrows():
            value = row['GeneSymbol|startPos|endPos|strand'].split('|')
            try:
                start = int(value[1])
            except:
                print(sampleID)
                print(i)
                print(pos)
                print(row)
                sys.exit(1)
            stop = int(value[2]) - start
            codon_pos = int(pos) - start
            codon_pos_percent = codon_pos / stop
            positions_in_ORF.append(codon_pos_percent)
            positions_in_triplet.append(pos%3 + 1)
            gene_symb = value[0]
            spy_gene_sequences[-1][gene_symb] = genome_seq[start:int(value[2])+1]

        
        ####the following doesn't make sense - need to fix the intent.
        
        # setup the plot
        ## using the triplet position
        fig, ax = plt.subplots()
        ax.set_xticks([1,2,3])
        ax.set_xticklabels(['1','2','3'])
        ax.set_title('$n_{ORFs}$ = '+str(len(positions_in_triplet)))
        n, bins, patches = plt.hist(positions_in_triplet, bins=3, rwidth=0.95,align='mid')
        plt.savefig(f"results/{sampleID[1:]}/{sampleID[1:]}.codon_position_{suffix}.png", dpi=300, format='png', bbox_inches='tight')
        
        ## Plotting a bar graph
        ## using the ORF position
        fig, ax = plt.subplots()
        ax.set_title('$n_{ORFs}$ = '+str(len(positions_in_triplet)))
        n, bins, patches = plt.hist(positions_in_ORF, bins=10, rwidth=0.95)
        plt.savefig(f"results/{sampleID[1:]}/{sampleID[1:]}.codon_position_ORF_{suffix}.png", dpi=300, format='png', bbox_inches='tight')
    
    return genes_ds, spy_gene_sequences

def run_cov_per_gene(genes_ds, df_sorted, sampleID):
    results = []
    for i in range(2):#pos, neg
        suffix = 'pos' if i == 0 else 'neg'
        
        ## Extract all genes and their names as lists
        genes_list = genes_ds[i]['GeneSymbol|startPos|endPos|strand'].to_list()
        alias_list = genes_ds[i]['GeneAlias|desc'].to_list()
        
        ## count occurance per gene
        ## and map gene to its alias in a separate list
        ## ... and get the sum of all reads per gene
        
        counts = []
        nr_reads = []
        aliases = []
        uniques = set(genes_list)
        
        for SPy in uniques:
            counts.append(genes_list.count(SPy))
            aliases.append(alias_list[genes_list.index(SPy)])
            spy_rows = df_sorted[i].loc[df_sorted[i]['GeneSymbol|startPos|endPos|strand'] == SPy]
            spy_sum = spy_rows[('Toxin', '+PNK')].sum()
            nr_reads.append(spy_sum)

        ## Now we make a new DateFrame from the three lists
        df = pd.DataFrame({'GeneID': list(uniques), 'Gene Name': aliases, "Site counts": counts, 'Sum reads': nr_reads})
        df_by_sites = df.sort_values(by=['Site counts','Sum reads'], ascending=False)
        df_by_sites.to_csv(f"results/{sampleID[1:]}/{sampleID[1:]}.df_by_sites.{suffix}.csv", index=True, header=True)
        # print(df_by_sites.head(25))

        results.append(df_by_sites)
    return results

def benjamini_hochberg(p_values):
    """
    Benjamini-Hochberg correction for multiple testing
    """
    n = len(p_values)
    sorted_indices = np.argsort(p_values)
    sorted_p_values = np.array(p_values)[sorted_indices]
    adjusted_p_values = np.empty(n, dtype=np.float64)

    cumulative_min = np.inf
    for i in range(n - 1, -1, -1):
        rank = i + 1
        adjusted_p = sorted_p_values[i] * n / rank
        cumulative_min = min(cumulative_min, adjusted_p)
        adjusted_p_values[sorted_indices[i]] = cumulative_min

    return adjusted_p_values

def run_pathway_analysis(sampleID, df_by_sites):
    print(F"Running pathway analysis for {sampleID[1:]}")
    for i in range(2):#pos, neg
        suffix = 'pos' if i == 0 else 'neg'
        
        geneIDs = list(df_by_sites[i]['GeneID'].apply(lambda x: x.split('|')[0].split(',')[-1].replace(' SPy', 'SPy')))
        filtered_geneIDs = [g for g in geneIDs if not g.startswith('SPY_RS')]
        
        ## R code - substituting with python equivalent
            # # Create the R script dynamically
            # r_script_content = f"""
            # library(clusterProfiler)
            # library(enrichplot)
            # library(ggplot2)

            # filtered_geneIDs <- c({', '.join([f'"{gene}"' for gene in filtered_geneIDs])})
            # print(length(filtered_geneIDs))
            # res <- enrichKEGG(filtered_geneIDs, organism = "spy", keyType = "kegg", pvalueCutoff = 0.5)
            # print(res$geneID)
            # print(res)
            # barplot(res, showCategory=10, title="Barplot of significantly enriched pathways")
            # barplot.savefig('barplot.png')
            # """

            # # Write the R script to a file
            # with open('run_analysis.R', 'w') as file:
            #     file.write(r_script_content)

            # # Execute the R script
            # result = subprocess.run(['Rscript', 'run_analysis.R'], capture_output=True, text=True)

            # # Print the output
            # print(result.stdout)
            # print(result.stderr)

        org_code = getattr(Defaults.Pairs, sampleID).kegg_org  # Replace with your organism code

        # Fetch pathway list for the organism
        if not os.path.exists('kegg'): os.makedirs('kegg')
        if not os.path.exists(f'kegg/kegg_pathways.{org_code}.txt'):
            print(F"Querying KEGG for pathway list", flush=True)
            pathway_list_url = f'http://rest.kegg.jp/list/pathway/{org_code}'
            pathway_list = requests.get(pathway_list_url).text
            # Parse the pathway list
            pathways = [line.split('\t') for line in pathway_list.strip().split('\n')]
            with open(f'kegg/kegg_pathways.{org_code}.txt', 'w') as outfile:
                for pathway in pathways:
                    outfile.write('\t'.join(pathway)+'\n')
        else:
            print(F"Reading local file kegg/kegg_pathways.{org_code}.txt", flush=True)
            with open(f'kegg/kegg_pathways.{org_code}.txt', 'r') as infile:
                pathways = [line.strip().split('\t') for line in infile]
        
        pathway_dict = {pathway[0].replace(f'path:{org_code}', ''): pathway[1] for pathway in pathways}

        # Fetch gene-pathway mappings
        if not os.path.exists(f'kegg/kegg_gene_pathways.{org_code}.txt'):
            with open(f'kegg/kegg_gene_pathways.{org_code}.txt', 'w') as outfile:
                print(F"Querying KEGG for gene-pathway mappings", flush=True)
                for pathway in pathway_dict.keys():
                    pathway_url = f'http://rest.kegg.jp/link/{org_code}/{pathway}'
                    gene_pathway_list = requests.get(pathway_url).text
                    for gp in gene_pathway_list.strip().split('\n'):
                        outfile.write(F"{pathway}\t{gp}\n")
        
        gene_pathway_dict = {}
        temp = {}
        with open(f'kegg/kegg_gene_pathways.{org_code}.txt', 'r') as infile:
            print(F"reading local file kegg/kegg_gene_pathways.{org_code}.txt", flush=True)
            for line in infile:
                line = line.strip().split('\t', 1)
                if not line[0] in temp: temp[line[0]] = []
                temp[line[0]].append(line[1])
        
        for pathway in pathway_dict.keys():
            gene_pathways = [line.split('\t') for line in temp[pathway]]
            for gene_pathway in gene_pathways:
                pathID = gene_pathway[0].replace(f'path:', '')
                gene = gene_pathway[1].replace(f'{org_code}:', '')
                if gene in gene_pathway_dict:
                    gene_pathway_dict[gene].append(pathID)
                else:
                    gene_pathway_dict[gene] = [pathID]

        # Convert the gene-pathway dictionary to a DataFrame
        gene_pathway_df = pd.DataFrame([(gene, pathway) for gene, pathways in gene_pathway_dict.items() for pathway in pathways], columns=['gene', 'pathway'])
        
        # Calculate the enrichment for each pathway
        pathway_enrichment = []
        M = len(gene_pathway_df['gene'].unique())  # Total number of unique genes
        N = len(filtered_geneIDs)  # Number of genes in the list

        
        for pathway, pathway_genes in gene_pathway_df.groupby('pathway'):
            pathway_genes_set = set(pathway_genes['gene'])
            k = len(pathway_genes_set)
            x = len(set(filtered_geneIDs) & pathway_genes_set)
            p_value = hypergeom.sf(x-1, M, k, N)
            
            ###following code is an alternative using fisher's exact test
            # # Create a contingency table
            # contingency_table = [[x, len(gene_list) - x],
            #                     [k - x, total_genes_in_pathways - k - (len(gene_list) - x)]]
            # # Perform Fisher's exact test
            # _, p_value = fisher_exact(contingency_table, alternative='greater')
            
            pathway_enrichment.append((pathway, pathway_dict[pathway], x, k, p_value))

        # Convert results to DataFrame
        enrichment_results = pd.DataFrame(pathway_enrichment, columns=['pathway', 'description', 'overlap', 'pathway_size', 'p_value'])

        # Apply Benjamini-Hochberg correction
        enrichment_results['adj_p_value'] = benjamini_hochberg(enrichment_results['p_value'].tolist())

        # Sort by p-value
        enrichment_results = enrichment_results.sort_values('adj_p_value')
        
        significant_pathways = enrichment_results[enrichment_results['adj_p_value'] < Defaults.pathway_enrich_adjust_p_threshold]
        # print(enrichment_results.head())
        
        # Plot the barplot of significantly enriched pathways
        plt.figure(figsize=(15, 8))
        sns.barplot(x='overlap', y='description', data=significant_pathways)
        plt.title('Barplot of Significantly Enriched Pathways')
        plt.xlabel('Number of Overlapping Genes')
        plt.ylabel('Pathway')
        plt.tight_layout()
        plt.savefig(f"results/{sampleID[1:]}/{sampleID[1:]}.pathway_analysis_{suffix}.png", dpi=300, format='png', bbox_inches='tight')
        print(F"Created figure:  results/{sampleID[1:]}/{sampleID[1:]}.pathway_analysis_{suffix}.png")

def pathway_enrich_calc_generic(gene_list, sampleID, suffix):
    org_code = getattr(Defaults.Pairs, sampleID).kegg_org  # Replace with your organism code

    # Fetch pathway list for the organism
    if not os.path.exists('kegg'): os.makedirs('kegg')
    if not os.path.exists(f'kegg/kegg_pathways.{org_code}.txt'):
        print(F"Querying KEGG for pathway list", flush=True)
        pathway_list_url = f'http://rest.kegg.jp/list/pathway/{org_code}'
        pathway_list = requests.get(pathway_list_url).text
        # Parse the pathway list
        pathways = [line.split('\t') for line in pathway_list.strip().split('\n')]
        with open(f'kegg/kegg_pathways.{org_code}.txt', 'w') as outfile:
            for pathway in pathways:
                outfile.write('\t'.join(pathway)+'\n')
    else:
        print(F"Reading local file kegg/kegg_pathways.{org_code}.txt", flush=True)
        with open(f'kegg/kegg_pathways.{org_code}.txt', 'r') as infile:
            pathways = [line.strip().split('\t') for line in infile]
    
    pathway_dict = {pathway[0].replace(f'path:{org_code}', ''): pathway[1] for pathway in pathways}

    # Fetch gene-pathway mappings
    if not os.path.exists(f'kegg/kegg_gene_pathways.{org_code}.txt'):
        with open(f'kegg/kegg_gene_pathways.{org_code}.txt', 'w') as outfile:
            print(F"Querying KEGG for gene-pathway mappings", flush=True)
            for pathway in pathway_dict.keys():
                pathway_url = f'http://rest.kegg.jp/link/{org_code}/{pathway}'
                gene_pathway_list = requests.get(pathway_url).text
                for gp in gene_pathway_list.strip().split('\n'):
                    outfile.write(F"{pathway}\t{gp}\n")
    
    gene_pathway_dict = {}
    temp = {}
    with open(f'kegg/kegg_gene_pathways.{org_code}.txt', 'r') as infile:
        print(F"reading local file kegg/kegg_gene_pathways.{org_code}.txt", flush=True)
        for line in infile:
            line = line.strip().split('\t', 1)
            if not line[0] in temp: temp[line[0]] = []
            temp[line[0]].append(line[1])
    
    for pathway in pathway_dict.keys():
        gene_pathways = [line.split('\t') for line in temp[pathway]]
        for gene_pathway in gene_pathways:
            pathID = gene_pathway[0].replace(f'path:', '')
            gene = gene_pathway[1].replace(f'{org_code}:', '')
            if gene in gene_pathway_dict:
                gene_pathway_dict[gene].append(pathID)
            else:
                gene_pathway_dict[gene] = [pathID]

    # Convert the gene-pathway dictionary to a DataFrame
    gene_pathway_df = pd.DataFrame([(gene, pathway) for gene, pathways in gene_pathway_dict.items() for pathway in pathways], columns=['gene', 'pathway'])
    
    # Calculate the enrichment for each pathway
    pathway_enrichment = []
    M = len(gene_pathway_df['gene'].unique())  # Total number of unique genes
    N = len(gene_list)  # Number of genes in the list
    
    print(gene_list[:5], len(gene_list))
    print(gene_pathway_df['gene'].unique()[:5], len(gene_pathway_df['gene'].unique()))
    print(len(set(gene_list) & set(gene_pathway_df['gene'].unique())))
    print(set(gene_list) - set(gene_pathway_df['gene'].unique()))
    sys.exit()
    
    for pathway, pathway_genes in gene_pathway_df.groupby('pathway'):
        pathway_genes_set = set(pathway_genes['gene'])
        k = len(pathway_genes_set)
        x = len(set(gene_list) & pathway_genes_set)
        p_value = hypergeom.sf(x-1, M, k, N)
        
        # ###following code is an alternative using fisher's exact test
        # # Create a contingency table
        # total_genes_in_pathways = len(pathway_genes['gene'].unique())
        # len_gene_list = len(set(gene_list)&set(pathway_genes['gene'].unique()))
        # contingency_table = [[x, len_gene_list - x],
        #                     [k - x, total_genes_in_pathways - k - (len_gene_list - x)]]
        # # Perform Fisher's exact test
        # _, p_value_fe = fisher_exact(contingency_table, alternative='greater')
        # print(F"{pathway} - {pathway_dict[pathway]} - {x} - {k} - {p_value} - {p_value_fe}")
        
        pathway_enrichment.append((pathway, pathway_dict[pathway], x, k, p_value))

    # Convert results to DataFrame
    enrichment_results = pd.DataFrame(pathway_enrichment, columns=['pathway', 'description', 'overlap', 'pathway_size', 'p_value'])

    # Apply Benjamini-Hochberg correction
    enrichment_results['adj_p_value'] = benjamini_hochberg(enrichment_results['p_value'].tolist())

    # Sort by p-value
    enrichment_results = enrichment_results.sort_values('adj_p_value')
        
    significant_pathways = enrichment_results[enrichment_results['adj_p_value'] < Defaults.pathway_enrich_adjust_p_threshold]
    
    # print(enrichment_results.head())
    # sys.exit()
    
    if len(significant_pathways) == 0:
        print(F"\nNo significant pathways found for {sampleID[1:]}\n")
    
    # Plot the barplot of significantly enriched pathways
    plt.figure(figsize=(15, 8))
    sns.barplot(x='overlap', y='description', data=significant_pathways)
    plt.title('Barplot of Significantly Enriched Pathways')
    plt.xlabel('Number of Overlapping Genes')
    plt.ylabel('Pathway')
    plt.tight_layout()
    outname = f"results/{sampleID[1:]}/{sampleID[1:]}.pathway_analysis.{suffix}.png"
    plt.savefig(outname, dpi=300, format='png', bbox_inches='tight')
    print(F"Created figure:  {outname}")

def download_kegg_pathways(organism='spy'):
    kegg_url = f"http://rest.kegg.jp/link/pathway/{organism}"
    response = requests.get(kegg_url)
    response.raise_for_status()  # Raise an error if the request was unsuccessful
    return response.text

def enrichment_analysis(gene_list, df_pathways):
    pathway_counts = df_pathways['Pathway'].value_counts()
    total_genes = len(df_pathways['Gene'].unique())
    results = []

    for pathway, count in pathway_counts.items():
        overlap = df_pathways[(df_pathways['Gene'].isin(gene_list)) & (df_pathways['Pathway'] == pathway)].shape[0]
        non_overlap = len(gene_list) - overlap
        pathway_genes = df_pathways[df_pathways['Pathway'] == pathway].shape[0]
        non_pathway_genes = total_genes - pathway_genes

        # Create contingency table
        contingency_table = np.array([[overlap, non_overlap],
                                      [pathway_genes - overlap, non_pathway_genes - non_overlap]])

        # Perform Fisher's exact test
        _, p_value = fisher_exact(contingency_table, alternative='greater')
        results.append((pathway, overlap, p_value))

    # Convert results to DataFrame
    df_results = pd.DataFrame(results, columns=['Pathway', 'Overlap', 'P-Value'])
    df_results = df_results.sort_values(by='P-Value').reset_index(drop=True)
    return df_results

def run_ident_cut_sites(genbank_file, sampleID, spy_gene_sequences):
    
    def get_best_align(in_seq):
        matrix = matlist.blosum62
        gap_open = -10
        gap_ext = -0.5
        best_score = {'name':'nan', 'score':0}
        score = 0

        counter = 0
        for name,source_seq in source_features.items():
            counter += 1
            aln_score = pairwise2.align.globalds(in_seq, str(source_seq),matrix, gap_open, gap_ext, score_only=True)
            if aln_score > score:
                score = aln_score
                best_score['name'] = name
                best_score['score'] = aln_score
                print(counter)
            if counter%100 == 0:
                print(counter)

        return best_score
    
    def align_worker(query_seq,source_gene,source_seq,q,threadID):
        matrix = matlist.blosum62
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -0.5
        aligner.extend_gap_score = -0.1
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0
        aln_score = pairwise2.align.globalms(source_seq,query_seq, 1,-.5, -.75, -.25, one_alignment_only=1, score_only=True)
        aln_score = aligner.score(source_seq,query_seq)
        q.put((source_gene,source_seq,aln_score,threadID))
        return(True)
    
    def calc_SD(liste):
        if len(liste)>1:
            avg=sum(liste)/len(liste)
            return (sum([(v-avg)**2 for v in liste])/(len(liste)-1))**.5
        else:
            return 'n/a'
        
    def find_homologs(source_seqs_orig, query_seq, cores, dev=False):
        source_seqs = copy.deepcopy(source_seqs_orig)
        q=Queue()
        thread_list=[]
        results=[]
        results2={}

        # thread counter
        count=0
        ticks = int(len(source_seqs)/100)
        print('|0%' + 92*'.' + '100%|')
        
        while len(source_seqs)>0 or len(thread_list)>0:
            while len(source_seqs)>0 and len(thread_list)<cores:
                count+=1
                if count%ticks == 0:
                    print('',end='|')
                source_gene,source_seq=source_seqs.popitem()
                thread_list.append((count,Process(target=align_worker,args=(query_seq,source_gene,source_seq,q,count))))
                thread_list[-1][1].deamon=True
                thread_list[-1][1].start()
            if not q.empty():
                source_gene,source_seq,score,threadID = q.get()
                results.append([score,source_gene,source_seq])
                #index == which thread was it
                index = [t[0] for t in thread_list].index(threadID)
                thread_list[index][1].join()
                thread_list.pop(index)

        results.sort(reverse=True)
        results2[results[0][1]] = {'score' : results[0][0], 'seq' : results[0][2]}
        results2[results[1][1]] = {'score' : results[0][0], 'seq' : results[0][2]}
        
        #to get a few negative cases for determing the proper threshold
        results=results[:30]
        avg=sum([s[0] for s in results])/len(results)
        SD=calc_SD([s[0] for s in results])
        while results[0][0]>=avg+2.5*SD:
            # Move results[0][1] over to results2
            # Then check the new results[0] (score)
            current_res = results.pop(0)
            results2[current_res[1]] = {'score' : current_res[0], 'seq': current_res[2]}
            
        ## Return list of everything that is at least
        ## 3*SD above the average (names only)
        if dev:
            return(results)
        else:
            return(results2)
    
    counter = 0
    source_features = {}
    source_feature_types = {}
    ref_sequence = ''
    for record in SeqIO.parse(genbank_file, 'genbank'):
        if ref_sequence == '':
            ref_sequence = record.seq.upper()
        for f in record.features[1:]:
            sequence = record.seq[f.location.start:f.location.end+1]
            input('check if code below is correct regarding strand <enter>')
            if f.location.strand != 1:
                sequence = reverse_complement(sequence)
        if 'locus_tag' not in f.qualifiers:
            locus_tag = 'no_tag_' + str(counter)
            counter += 1
        else:
            locus_tag = f.qualifiers["locus_tag"][0]
        source_features[locus_tag] = str(sequence)
        source_feature_types[locus_tag] = f.type
    
    ## Iterate through the identified sequences
    homologs = [{}, {}]#pos, neg
    
    for i in range(2):#pos, neg
        suffix = 'pos' if i == 0 else 'neg'
        counter = 0
        for ID, seq in spy_gene_sequences[i].items():
            counter += 1
        # print('\n'+str(counter),"ID:", ID)
        
        res = find_homologs(source_features, ref_sequence, cores=Defaults.cores)
        if len(res) > 0:
            homologs[i][ID] = res
        # print(res.keys())
        
        ##saving results to json file
        with open(F"results/{sampleID[1:]}/{sampleID[1:]}.homologs_{suffix}.json", 'w') as outfile:
            json.dump(homologs, outfile)

    ## get aliases from KEGG
    url = 'http://rest.kegg.jp/get/lmo:lmo2558'
    r = requests.get(url)
    lines = r.text.split('\n')

    infos = {'exp_ID':[],
            'exp_NAME' : [],
            'exp_ORTH' : [],
            'source_ID' : [],
            'source_NAME' : [],
            'source_ORTH' : []}


    ## Use the known info for the experiement organism
    df_homol = pd.DataFrame(columns=["GeneID",'Gene Name','Source GeneID'])
    for exp_id in homologs:
        df_homol = df_homol.append(df_by_sites.loc[df_by_sites["GeneID"].str.contains(exp_id)].reset_index()[["GeneID",'Gene Name']])
        df_homol["GeneID"] = df_homol["GeneID"].apply(lambda x : x.split('|')[0])
        df_homol["Source GeneID"] = df_homol["GeneID"].apply(lambda x : list(homologs[x.split('|')[0]].keys()))

    return source_features

def run_motif_search(genome_seq, sampleID):
    motif_re = ambiguity_to_regex(Defaults.motif)
    cutsites = []
    ## Whole genome
    for match in re.finditer(motif_re, str(genome_seq)):
        cutsite = match.span()[0]+2
        cutsites.append(cutsite)
    print("Total number of cutsites in the genome using the found motif:",len(cutsites))

    cutsites_f = 0
    cutsites_t = 0
    cutsites_cds = 0

    ## Annotated features
    cutsites_dict = {}
    cutsites_cds_dict = {}
    for locus_tag,seq in source_features.items():
        iterator = re.finditer(regex,str(seq))
        cutsites_this = 0
        for match in iterator:
            #cutsite = match.span()[0]+2
            cutsites_f += 1
            if re.match('no_tag',locus_tag) is None:
                cutsites_t += 1
                cutsites_this += 1
            if source_feature_types[locus_tag] == 'CDS':
                cutsites_cds += 1
        if source_feature_types[locus_tag] == 'CDS' and cutsites_this > 0:       
            cutsites_cds_dict[locus_tag] = [cutsites_this]
        cutsites_dict[locus_tag] = [cutsites_this]
    print("Total number of cutsites in annotated features using the found motif:",cutsites_f)
    print("Total number of cutsites in named features using the found motif:",cutsites_t)
    print("Total number of cutsites in CDS using the found motif:",cutsites_cds)
    print()

    df_cutsites = pd.DataFrame.from_dict(cutsites_dict, orient='index').sort_values(by=0, ascending=False)
    #print(df_cutsites.head())

    feature_w_cutsites = df_cutsites[0].loc[df_cutsites[0] > 0].to_list()
    cds_w_cutsites = len(cutsites_cds_dict)
    
    proportion = len(feature_w_cutsites)/len(source_features)
    print('Proportion of all features:',round(proportion*100),'%')
    proportion = cds_w_cutsites/len([t for t in source_feature_types if source_feature_types[t] == 'CDS'])
    print('Proportion of CDS:',round(proportion*100),'%')
    proportion = len([cs for cs in cutsites_cds_dict if cutsites_cds_dict[cs][0] >= 3])/len([t for t in source_feature_types if source_feature_types[t] == 'CDS'])
    print('Proportion of CDS with min. 3 cutsites:',round(proportion*100),'%')

def select_genes_w_motifs(genbank_file, sampleID):
    '''assuming that each gb file contains only one genome'''
    motif_re = ambiguity_to_regex(Defaults.motif)
    
    ref_sequence = ''
    counter = 0
    dataset = []
    
    for record in SeqIO.parse(genbank_file, 'genbank'):
        # if ref_sequence == '':
        #     ref_sequence = record.seq.upper()
        for feature in record.features:
            if feature.type != 'source':
                sequence = record.seq[feature.location.start:feature.location.end+1]
                if feature.location.strand == -1:
                    sequence = sequence.reverse_complement()
                if 'old_locus_tag' in feature.qualifiers:
                    locus_tag = feature.qualifiers['old_locus_tag'][0]
                elif 'locus_tag' in feature.qualifiers:
                    locus_tag = feature.qualifiers['locus_tag'][0]
                else:
                    locus_tag = 'no_tag_' + str(counter)
                    counter += 1
                cut_count = len(re.findall(motif_re, str(sequence)))
                dataset.append((feature.type, locus_tag, cut_count, feature.location.start, feature.location.end+1))

    total_CDS = len([d for d in dataset if d[0] == 'CDS'])
    min_1_cutsites = [d for d in dataset if d[2] > 0 and d[0] == 'CDS']
    min_3_cutsites = [d for d in dataset if d[2] >= 3 and d[0] == 'CDS']

    print(F"there are {sum([d[2] for d in dataset])} cutsites in total")
    print(F"there are {len(min_1_cutsites)} CDS with at least one cutsites {len(min_1_cutsites)/total_CDS*100:.1f}%")
    print(F"there are {len(min_3_cutsites)} CDS with at least three cutsites {len(min_3_cutsites)/total_CDS*100:.1f}%")

    return [d for d in dataset if d[2] >= Defaults.min_motif_per_CDS and d[0] == 'CDS']

def reorder_sam(samfile, fastq, output_sam):
    read_order = []
    print(F"Reading {fastq}", flush=True)
    with open(fastq, 'r') as infile:
        for e,line in enumerate(infile):
            if e%4 == 0:
                read_order.append(line.strip().split()[0][1:])
                if e % (4*1000000) == 0:
                    print(F"Read {e/4:,} seq IDs", end='\r', flush=True)
    print()
    
    sam_header = []
    sam_content = {}
    print(F"Reading {samfile}", flush=True)
    with open(samfile, 'r') as infile:
        for e,line in enumerate(infile):
            if line.startswith('@'):
                sam_header.append(line.strip())
            else:
                read_id = line.split('\t')[0]
                if not read_id in sam_content:
                    sam_content[read_id] = [line.strip()]
                else:
                    sam_content[read_id].append(line.strip())
            if e%1000000 == 0:
                print(F"Read {e:,} lines", end='\r', flush=True)
    print()
    
    print(F"Writing {output_sam}", flush=True)
    count = 0
    with open(output_sam, 'w') as outfile:
        for line in sam_header:
            outfile.write(line+'\n')
        for read_id in read_order:
            if read_id in sam_content:
                for line in sorted(sam_content[read_id]):
                    outfile.write(line+'\n')
                    count += 1
                    if count % 1000000 == 0:
                        print(F"Wrote {count:,} lines", end='\r', flush=True)
    print()
    return 1

def create_new_sample_pair(sampleID, subrefID, subsetID, refseq, sample_datasets):
    new_sampleID = F"{sampleID}_{subsetID}"

    if not hasattr(Defaults.Pairs, new_sampleID):
        setattr(Defaults.Pairs, new_sampleID, sample_data())
        for key, val in sample_datasets[sampleID].items():
            setattr(getattr(Defaults.Pairs, new_sampleID), key, val)
        print(F"Initialized {new_sampleID}")
    
    
    print(F"Preparing reference data for genome subset {subsetID} of {sampleID} with subrefID {subrefID}")
    
    #ref gff
    ref_set = set()
    setattr(getattr(Defaults.Pairs, new_sampleID), 'name_of_gene_information_file', F"references/{new_sampleID}.{subrefID}.gene_information.csv")
    # if not os.path.exists(F"references/{new_sampleID}.{subrefID}.gene_information.csv"):
    print(F"Creating gene information file for {new_sampleID}")
    with open(F"references/{new_sampleID}.{subrefID}.gene_information.csv", 'w') as outfile:
        with open(getattr(Defaults.Pairs, sampleID).name_of_gene_information_file, 'r') as infile:
            for line in infile:
                if line.startswith('#'):
                    outfile.write(line)
                else:
                    ref_set.add(line.split()[0])
                if line.startswith(subrefID):
                    outfile.write(line)
    print(F"References found in {sampleID}: {', '.join(ref_set)}")
    #validating separation
    ref_set = set()
    with open(F"references/{new_sampleID}.{subrefID}.gene_information.csv", 'r') as infile:
        for line in infile:
            if len(line)>2 and line[0] != '#':
                ref_set.add(line.split()[0])
    print(F"References found in {new_sampleID}: {', '.join(ref_set)}")
    
    if getattr(Defaults.Pairs, sampleID).name_of_gene_information_file == getattr(Defaults.Pairs, new_sampleID).name_of_gene_information_file:
        print(F"Gene information file for {new_sampleID} is the same as for {sampleID}")
        print(getattr(Defaults.Pairs, sampleID).name_of_gene_information_file)
        sys.exit(1)
    
    #ref seq
    setattr(getattr(Defaults.Pairs, new_sampleID), 'name_of_genome_file', F"references/{new_sampleID}.{subrefID}.fasta")
    # if not os.path.exists(F"references/{new_sampleID}.{subrefID}.fasta"):
    print(F"Creating fasta file for {new_sampleID}")
    with open(F"references/{new_sampleID}.{subrefID}.fasta", 'w') as outfile:
        outfile.write(F">{subrefID}\n{refseq}\n")
        
    return new_sampleID

def get_number_of_reads(fastq):
    count = 0
    plus = False
    with open(fastq, 'r') as infile:
        for line in infile:
            if line.startswith('@') and not plus:
                count += 1
                plus = True
            elif line.startswith('+'):
                plus = False
        return count

def init_sample_datasets(sample_datasets):
    for sampleID in sample_datasets:
        if not hasattr(Defaults.Pairs, sampleID):
            setattr(Defaults.Pairs, sampleID, sample_data())
            for key, val in sample_datasets[sampleID].items():
                setattr(getattr(Defaults.Pairs, sampleID), key, val)
            print(F"Initialized {sampleID}")

def extract_read_count_from_sam(samfile):
    count = []
    with open(samfile, 'r') as infile:
        for line in infile:
            if not line.startswith('@'):
                if len(line)>2:
                    count.append(line.split()[0])
    return len(set(count))

def filter_sam_file(sampleID, refID):
    content = []
    counts = [0, 0]
    with open(getattr(Defaults.Pairs, sampleID).name_of_sam_file, 'r') as infile:
        for line in infile:
            counts[0] += 1
            if line.startswith('@'):
                content.append(line)
            else:
                if line.split()[2] == refID:
                    content.append(line)
    with open(getattr(Defaults.Pairs, sampleID).name_of_sam_file, 'w') as outfile:
        for line in content:
            counts[1] += 1
            outfile.write(line)
    
    print(F"Filtered {sampleID} sam file from {counts[0]:,} to {counts[1]:,} lines")
    
    return True

def main(sample_datasets):
    
    #initializing sample datasets
    init_sample_datasets(sample_datasets)
        
    with open(Defaults.log_file, 'w') as logfile:
        # #calculate reference indices if needed
        logfile.write('Calculating reference indices\n')
        create_reference_index()
        
        #running perl script for emote analysis
        for sampleID in dir(Defaults.Pairs):
            if not sampleID.startswith('__') and sampleID == 'PSF370_A1_phage_30min_only_A1ref':
                # if not os.path.exists('results/'+sampleID[1:]):
                    logfile.write(F"Running emote analysis for {sampleID}\n")
                    print(f"Sample ID: {sampleID}")
                    
                    logfile.write(F"Reading ref fasta for {sampleID}\n")
                    genome_seqs = read_fasta(getattr(Defaults.Pairs, sampleID).name_of_genome_file)
                    logfile.write(F"Read ref fasta for {sampleID}\n")

                    Defaults.Definitions = read_replicate_definition(Defaults.Definitions_file, sampleID)

                    #handling of multiple reference genomes
                    subsetID = 0
                    sampleID_iteration = sampleID
                    for refID, refseq in genome_seqs:
                        subsetID += 1
                        if len(genome_seqs) > 1: 
                            print(f'Splitting up sample into subsets, processing #{subsetID}')
                            #creating new sample pair object
                            sampleID_iteration = create_new_sample_pair(sampleID, refID, subsetID, refseq, sample_datasets)
                            
                            create_reference_index(target=getattr(Defaults.Pairs, sampleID_iteration).name_of_genome_file)
                        
                        run_emote_step1(sampleID_iteration)
                        logfile.write(F"Ran emote step 1 for {sampleID_iteration}\n")
                        
                        #check and report number of reads
                        pre_count = get_number_of_reads(getattr(Defaults.Pairs, sampleID_iteration).name_of_the_fastq_file)
                        post_count = get_number_of_reads(getattr(Defaults.Pairs, sampleID_iteration).name_of_the_output_file+'Emote.fastq')
                        logfile.write(F"Number of reads for {sampleID_iteration} before emote step 1: {pre_count:,}\n")
                        logfile.write(F"Number of reads for {sampleID_iteration} after emote step 1: {post_count:,}\n")
                    
                        if post_count == 0:
                            print('No read remains after emote step 1')
                            logfile.write(F"Skipping {sampleID_iteration} due to 0 reads after emote step 1\n")
                            
                            print(sampleID, getattr(Defaults.Pairs, sampleID).name_of_the_fastq_file)
                            for i in range(1, subsetID+1):
                                print(F"{sampleID}_{i}", getattr(Defaults.Pairs, F"{sampleID}_{i}").name_of_the_fastq_file)
                            sys.exit('Error')
                        else:
                    
                            run_bowtie_mapping(sampleID_iteration, sampleID)
                            logfile.write(F"Ran bowtie mapping for {sampleID_iteration}\n")
                            
                            #filter sam file accordingly
                            filter_sam_file(sampleID_iteration, refID)
                            
                            run_emote_step2(sampleID_iteration)
                            logfile.write(F"Processed sam file for {sampleID_iteration} (emote step 2)\n")
                            
                            read_count = extract_read_count_from_sam(getattr(Defaults.Pairs, sampleID_iteration).name_of_sam_file)
                            logfile.write(F"Dataset {sampleID_iteration} has {read_count:,} reads in mapping file\n")
                            
                            df_sorted = run_follow_up_analysis(sampleID_iteration)#pos, neg
                            logfile.write(F"Ran follow up analysis for {sampleID_iteration}\n")
                            try:
                                run_motif_scan(refseq, df_sorted, sampleID_iteration)
                                logfile.write(F"Ran motif scan for {sampleID_iteration}\n")
                            except:
                                logfile.write(F"Failed to run motif scan for {sampleID_iteration}\n")
                                print('Failed to run motif scan (completely)')
                            
                            try:
                                genes_ds, spy_gene_sequences = run_codon_nucleotide_position(df_sorted, genome_seq, sampleID_iteration)
                                logfile.write(F"Ran codon nucleotide position for {sampleID_iteration}\n")
                                df_by_sites = run_cov_per_gene(genes_ds, df_sorted, sampleID_iteration)
                                logfile.write(F"Ran coverage per gene for {sampleID_iteration}\n")
                                run_pathway_analysis(sampleID_iteration, df_by_sites)
                                logfile.write(F"Ran pathway analysis for {sampleID_iteration}\n")
                            except:
                                logfile.write(F"Failed to run codon nuc position, cov per gene, or pathway analysis for {sampleID_iteration}\n")
                                print('Failed to run pathway analysis')
                            """
                            source_features = run_ident_cut_sites(getattr(Defaults.Pairs, sampleID).genbank_file, sampleID, spy_gene_sequences)
                            run_motif_search(genome_seq, source_features, sampleID)
                            # CDS_w_cutsites = select_genes_w_motifs(getattr(Defaults.Pairs, sampleID).genbank_file, sampleID) #--> not being used at all anymore
                            # pathway_enrich_calc_generic([d[1] for d in CDS_w_cutsites], sampleID, 'CDS_w_3_cutsites')
                            """
                            # sys.exit('checkpoint')
    
    
    
    
    
#

main(sample_datasets)

print('Done')