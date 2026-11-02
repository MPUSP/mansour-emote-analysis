#!/usr/bin/env python3

import os
import sys
import time

from argparse import ArgumentParser
parser=ArgumentParser()
parser.add_argument('--fq', dest='fastq', default='', help='read fastq file')
parser.add_argument('-s', dest='samfile', default='', help='sam alignment file')
parser.add_argument('-r', dest='ref', default='', help='reference genome fasta file')
parser.add_argument('-a', dest='annotation', default='', help='reference annotation file (gff3)')
parser.add_argument('-b', dest='barcode', default='', help='emote barcode file')
parser.add_argument('-o', dest='prefix', default='./', help='output prefix')
args=parser.parse_args()#global variable

def check_input():
    if args.fastq == '':
        print('Please provide a fastq file')
        sys.exit(1)
    elif not os.path.exists(args.fastq):
        print(f'FastQ file "{args.fastq}" does not exist')
        sys.exit(1)
        
    if args.samfile == '':
        print('Please provide a sam file')
        sys.exit(1)
    elif not os.path.exists(args.samfile):
        print(f'SAM file "{args.samfile}" does not exist')
        sys.exit(1)
    
    if args.ref == '':
        print('Please provide a reference genome file')
        sys.exit(1)
    elif not os.path.exists(args.ref):
        print(f'Reference genome file "{args.ref}" does not exist')
        sys.exit(1)
        
    if args.barcode == '':
        print('Please provide a barcode file')
        sys.exit(1)
    elif not os.path.exists(args.barcode):
        print(f'Barcode file "{args.barcode}" does not exist')
        sys.exit(1)
        
    if args.prefix == '': args.prefix = './'
    elif os.path.exists(args.prefix) and os.path.isdir(args.prefix) and not args.prefix.endswith('/'):
        args.prefix += '/'
    
    return True
        
def read_barcodes(barcode_file):
    dataset = {}
    with open(barcode_file, 'r') as f:
        header = False
        for line in f:
            if not header:
                header = True
                continue
            else:
                fields = line.strip().split('\t')#sampleID, BC seq
                dataset[fields[0]] = [fields[1], 0, 0]#BC seq, mapped reads per emote barcode, total number of reads
    return dataset

def read_fastq(fastq_file, emote_barcodes):#just to get the link between read ID and sample ID (barcode), also count total reads per bc
    dataset = {}
    count = 0
    with open(fastq_file, 'r') as f:
        entry = []
        for line in f:
            entry.append(line.strip())
            if len(entry) == 4:
                bc =  entry[2][1:].split(' ')[0]
                quantif = entry[2][1:].split(' ')[1]
                dataset[entry[0][1:].split(' ')[0]] = (bc, quantif)
                entry = []
                if count % 1000000 == 0:
                    print(F"Processed {count:,} fastq reads", end='\r')
                count += 1
                emote_barcodes[bc][2] += 1
    print()
    return dataset, emote_barcodes

def print_total_runtime(start_time, intermediate=False):
    end_time = time.time()
    total_time = end_time - start_time
    hours = total_time // 3600
    minutes, seconds = divmod(total_time % 3600, 60)
    if intermediate:
        print(F'Intermediate total runtime: {int(hours)}:{int(minutes)}:{seconds:.1f} H:M:S.')
    else:
        print(F'Total runtime: {int(hours)}:{int(minutes)}:{seconds:.1f} H:M:S.')
    return True

def read_fasta(fasta_file, choice=0):#choice 0-based, to select from multiple sequences fasta file
    dataset = []
    with open(fasta_file, 'r') as infile:
        title = ''
        seq = ''
        for line in infile:
            if line[0] == '>' and title == '':
                title = line.strip()[1:]
            elif line[0] == '>' and title != '':
                dataset.append([title, seq])
                title = line.strip()[1:]
                seq = ''
                if len(dataset) == choice + 1:
                    break
            else:
                seq += line.strip().upper()
        dataset.append([title, seq])
    
    return dataset[choice]

def read_annotation(anno_file):
    dataset = {}
    anno_lookup = {}
    with open(anno_file) as f:
        for line in f:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                if line[1] == 'feature' and line[2] != 'source':
                    key = (int(line[3]), int(line[4]), line[6])
                    if not line[0] in anno_lookup:
                        anno_lookup[line[0]] = {'+': {}, '-': {}}
                        dataset[line[0]] = {}

                    if not key in dataset[line[0]]:
                        dataset[line[0]][key] = {}
                        
                    line[8] = dict(e.split('=') for e in line[8].split(';'))
                    dataset[line[0]][key].update(line[8])
                    
                    for p in range(int(line[3]), int(line[4]) + 1):
                        if not p in anno_lookup[line[0]][line[6]]:
                            anno_lookup[line[0]][line[6]][p] = [key]
                        else:
                            anno_lookup[line[0]][line[6]][p].append(key)

    return dataset, anno_lookup
    
def convert_bitflag(number):
    bitflag = bin(int(number))[2:][::-1]
    bitflag = bitflag + '0'*(12-len(bitflag))
    return bitflag

def reverse_complement(string):
  table=bytes.maketrans(b'ACGTNWSMKRYBDHV.acgtnwsmkrybdhv',b'TGCANWSKMYRVHDBNtgcanwskmyrvhdb')
  string=string.encode('UTF-8').translate(table)[::-1]
  return(string.decode('unicode_escape'))
    
def processing_sam_file(sam_file, readID_bc, emote_barcodes):
    dataset = {}
    # quickly identifying read duplicates
    duplicates = {}
    with open(sam_file, 'r') as f:
        for line_count, line in enumerate(f):
            if line_count % 1000000 == 0:
                print(F"Processed {line_count:,} SAM lines for duplicated reads identification", end='\r')
            if not line.startswith('@'):
                line = line.strip().split('\t')
                bitflag = convert_bitflag(line[1])
                if not bitflag[2] == '1':#hence read was mapped
                    readID = line[0]
                    bc, quantif = readID_bc[readID]
                    if not bc in duplicates: duplicates[bc] = {}
                    if readID in duplicates[bc]:
                        duplicates[bc][line[0]].add(line_count)
                    else:
                        duplicates[bc][line[0]] = set([line_count])
    print(F"Processed {line_count:,} SAM lines for duplicated reads identification")
    
    dup_indeces = {}
    for bc in duplicates:
        if not bc in dup_indeces: dup_indeces[bc] = set()
        for readID in duplicates[bc]:
            if len(duplicates[bc][readID]) > 1:
                dup_indeces[bc] |= duplicates[bc][readID]
            emote_barcodes[bc][1] += 1#counting mapped reads
    
    for bc in emote_barcodes:
        if not bc in dup_indeces: dup_indeces[bc] = set()
    
    del duplicates
    
    count_mapped = {}



    with open(sam_file, 'r') as f:
        for line_count, line in enumerate(f):
            if line_count % 1000000 == 0:
                print(F"Processed {line_count:,} SAM lines for ends identification", end='\r')
            if not line.startswith('@'):
                line = line.strip().split('\t')
                readID = line[0]
                bitflag = convert_bitflag(line[1])
                if not bitflag[2] == '1':#hence read was mapped
                    bc, quantif = readID_bc[readID]
                    strand = 'plus' if bitflag[4] == '0' else 'minus'
                    if strand == 'plus': pos = int(line[3])
                    else: pos = int(line[3]) + len(line[9]) - 1#this method does not allow for InDels
                    if not pos in dataset:
                        dataset[pos] = {}
                    if not strand in dataset[pos]:
                        dataset[pos][strand] = {}
                    if not bc in dataset[pos][strand]:
                        dataset[pos][strand][bc] = [[], []]#unambiguous, ambiguous
                    if line_count in dup_indeces[bc]:
                        dataset[pos][strand][bc][1].append(quantif)
                    else:
                        dataset[pos][strand][bc][0].append(quantif)
                     
                    if not bc in count_mapped:
                        count_mapped[bc] = set()
                    count_mapped[bc].add(readID)   

    for bc in count_mapped:
        emote_barcodes[bc][1] = len(count_mapped[bc])
                        # if strand == 'plus': prim5_base = ref_seq[pos - 1]#that can be calculated in the end, waste of time here
                        # else: prim5_base = reverse_complement(ref_seq[pos + 1])
                        # prim5_base = prim5_base.replace('T', 'U')
    print(F"Processed {line_count:,} SAM lines for ends identification")

    return dataset, emote_barcodes

def get_nucleotide(pos, strand, ref_seq):
    if strand == 'plus': 
        base = ref_seq[pos - 1]
    else: 
        base = reverse_complement(ref_seq[pos + 1])
    if base == 'T': base = 'U'
    return base

def get_gene_info(pos, strand, refID, ref_annotation, anno_lookup):
    gene_info = []
    strand_symbol = '+' if strand == 'plus' else '-'
    if pos in anno_lookup[refID][strand_symbol]:
        for key in anno_lookup[refID][strand_symbol][pos]:
            c1 = ""
            if 'locus_tag' in ref_annotation[refID][key]:
                c1 = ref_annotation[refID][key]['locus_tag']+'|'
            else:
                c1 = 'n/a|'
            c1 += f"{key[0]}|{key[1]}|{'plus' if key[2] == '+' else 'minus'}"
            
            if 'gene' in ref_annotation[refID][key]:
                c2 = ref_annotation[refID][key]['gene']+'|'
            elif 'old_locus_tag' in ref_annotation[refID][key]:
                c2 = ref_annotation[refID][key]['old_locus_tag']+'|'
            elif 'db_xref' in ref_annotation[refID][key]:
                c2 = ref_annotation[refID][key]['db_xref']+'|'
            else:
                c2 = "n/a|"
            
            if 'GO_function' in ref_annotation[refID][key]:
                c2 += ref_annotation[refID][key]['GO_function']
            elif 'product' in ref_annotation[refID][key]:
                c2 += ref_annotation[refID][key]['product']
            else:
                c2 += "n/a"
                
            gene_info.append([c1, c2])
    return gene_info

def create_output_tables(emote_barcodes, results, refID, ref_seq, ref_annotation, anno_lookup):
    outpath = args.prefix
    if not outpath[-1] == '/' and not outpath[-1] == '.': outpath += '.'
    
    barcode_report = []
    
    with open(F"{outpath}EmoteBarcodesReport.csv", 'w') as f:
        f.write(f'Emote_Barcode\tTotal_emote_barcode_reads\tMapped_reads\tUnmapped_reads\tPercentage_mapped_reads\tPercentage_unmapped_reads\n')
        totals = [0, 0, 0]
        for bc in emote_barcodes:
            total_reads = emote_barcodes[bc][2]
            mapped_reads = emote_barcodes[bc][1]
            unmapped_reads = total_reads - mapped_reads
            totals[0] += total_reads
            totals[1] += mapped_reads
            totals[2] += unmapped_reads
            perc_mapped = mapped_reads / total_reads * 100
            perc_unmapped = 100 - perc_mapped
            f.write(F"{bc}\t{total_reads}\t{mapped_reads}\t{unmapped_reads}\t{perc_mapped:.2f}%\t{perc_unmapped:.2f}%\n")
        f.write(F"\nTotal\t{totals[0]}\t{totals[1]}\t{totals[2]}\n")
    print(F"Emote barcode report written to {outpath}EmoteBarcodesReport.csv")
    
    barcode_order = sorted(list(emote_barcodes.keys()))
    
    for strand in ['plus', 'minus']:
        for status in ['amb', 'unamb']:
            status_index = 1 if status == 'amb' else 0
            strand_name = 'Pos' if strand == 'plus' else 'Neg'
            with open(F"{outpath}{status.capitalize()}{strand_name}Table.csv", 'w') as f:
                f.write("Genomic_pos\tStrand\t5'_base\t")
                for bc in barcode_order:
                    f.write(F"{bc}\t")
                for bc in barcode_order:
                    f.write(F"{bc}_Q\t")
                f.write('GeneSymbol|startPos|endPos|strand\tGeneAlias|desc\n')
                for pos in sorted(results):
                    if strand in results[pos]:
                        line = [str(pos), strand, get_nucleotide(pos, strand, ref_seq)]
                        line_ctr = 0
                        for bc in barcode_order:
                            if bc in results[pos][strand]:
                                line.append(str(len(results[pos][strand][bc][status_index])))
                                line_ctr += len(results[pos][strand][bc][status_index])
                            else:
                                line.append('0')
                        for bc in barcode_order:
                            if bc in results[pos][strand]:
                                line.append(str(len(set(results[pos][strand][bc][status_index]))))
                            else:
                                line.append('0')
                        
                        if line_ctr > 0:
                            gene_info = get_gene_info(pos, strand, refID, ref_annotation, anno_lookup)
                            if len(gene_info) == 0:
                                f.write('\t'.join(line+['',''])+'\n')
                            else:
                                for gi in gene_info:
                                    f.write('\t'.join(line+gi)+'\n')
            print(F"Output table {outpath}{status.capitalize()}{strand.capitalize()}Table.csv written.")
        
    return True

def main():
    start_time_total = time.time()
    
    print('Checking input...')
    check_input()
    
    print(F"Reading reference sequence file {args.ref}")
    refID, ref_seq = read_fasta(args.ref, 0)#to select the first sequence in the fasta file
    refID = refID.split(' ')[0]
    
    print(F"Reading barcode file {args.barcode}")
    emote_barcodes = read_barcodes(args.barcode)
    
    print(F"Reading fastq file {args.fastq}")
    readID_bc, emote_barcodes = read_fastq(args.fastq, emote_barcodes)
    
    print(F"Reading reference annotation file {args.annotation}")
    ref_annotation, anno_lookup = read_annotation(args.annotation)
    
    print(F"Setting up additional variables")
    gene_info_col_tab = {"plus": {0:1, 1:1}, "minus": {0:1, 1:1}}# initialize the variables of the gene info column headers in the 4 tables.
    total_reads = sum([emote_barcodes[bc][2] for bc in emote_barcodes])
    
    print_total_runtime(start_time_total, intermediate=True)
    
    print(F"Processing SAM file {args.samfile}")
    results, emote_barcodes = processing_sam_file(args.samfile, readID_bc, emote_barcodes)
    
    print(F"Creating output tables")
    create_output_tables(emote_barcodes, results, refID, ref_seq, ref_annotation, anno_lookup)
    
    print_total_runtime(start_time_total)
    
main()

