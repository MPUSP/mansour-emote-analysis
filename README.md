# Analysis code for "A non-canonical type III toxin–antitoxin system encoded in *Streptococcus pyogenes* prophage promotes phage resistance"

Code accompanying the manuscript (under revision at PNAS; DOI to be added upon publication):

> Moïse Mansour, Chu Wang, Karin Hahnke, Fabián A. Cornejo, Michael Jahn, Kathirvel Alagesan, Sebastian Rämisch, Dior Beerens, Rina Ahmed-Begrich, Katja Schmidt, Florian Kondrot, Knut Finstermeier, Emmanuelle Charpentier.
> *A non-canonical type III toxin–antitoxin system encoded in Streptococcus pyogenes prophage promotes phage resistance.*

The repository contains the scripts used to analyse EMOTE-seq (Exact Mapping Of Transcriptome Ends) data and to generate the figures of the publication. Raw sequencing data are available at the European Nucleotide Archive under study accession **PRJEB109274**.

## Contents

| File | Description |
|------|-------------|
| `analysis.py` | Main analysis: runs the EMOTE analysis and creates the figures required for the manuscript. |
| `emote_step_II.py` | Python replacement for `emoteStepII.pl` of EMOTE-conv, with a fix for the parallel-processing issue of the original Perl code. |
| `LICENSE` | MIT license. |

## Inputs

- Sequencing reads (FASTQ)
- Reference genome sequences (FASTA)
- Reference annotations (GFF)

The EMOTE-seq reads (single-end, Illumina MiSeq) are deposited at ENA under **PRJEB109274**:

| Library | Sample accession |
|---------|------------------|
| Nemote-spy0938-Qz-Ecoli_S1 | ERS29556842 |
| Nemote-spy0938-Qz_S1 | ERS29556843 |
| Phage-30min-gg_S1 | ERS29556844 |
| EMOTE-phage-mut-aI_S1 | ERS29556845 |

## `emote_step_II.py`

Python version of `emoteStepII.pl` from **EMOTE-conv**, which converts EMOTE data into lists of quantified genomic 5'-end positions with associated genomic information. It is based on:

> Redder P, Yasrebi H (2015). EMOTE-conv: A Computational Pipeline to Convert Exact Mapping of Transcriptome Ends (EMOTE) Data to the Lists of Quantified Genomic Positions Correlated to Related Genomic Information. *J Appl Bioinforma Comput Biol* 4:118. doi:[10.4172/2329-9533.1000118](http://dx.doi.org/10.4172/2329-9533.1000118)

The original Perl code (<http://sourceforge.net/projects/emoteconv>) had a parallel-processing issue, which is fixed here. Please also cite the original EMOTE-conv publication when using this script.

```bash
python emote_step_II.py --fq reads.fastq -s alignment.sam -r genome.fasta \
    -a annotation.gff3 -b barcodes.txt -o output_prefix
```

| Argument | Description |
|----------|-------------|
| `--fq` | Read FASTQ file |
| `-s` | SAM alignment file |
| `-r` | Reference genome FASTA file |
| `-a` | Reference annotation file (GFF3) |
| `-b` | EMOTE barcode file |
| `-o` | Output prefix (default: `./`) |

## `analysis.py`

The full analysis. All input and output paths are hardcoded and refer to the authors' file system, so adjust the paths at the top of the script to your setup. There is no command-line interface; the script is simply executed:

```bash
python analysis.py
```

## Requirements

There is no installation routine. The code was run with Python 3.9 and the following packages, installed with pip:

`weblogo`, `logomaker`, `biopython` (uses `Bio.pairwise2`, which is not available in recent Biopython releases), `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`, `requests`

`analysis.py` additionally calls the following external programs via `subprocess`:

- [Bowtie](https://bowtie-bio.sourceforge.net/) 1.1.2 (`bowtie-build` and `bowtie`) for building the reference index and mapping reads
- Perl and the original [EMOTE-conv](http://sourceforge.net/projects/emoteconv) scripts, expected in a `scripts/` directory next to `analysis.py`: `emoteStepI.pl` with its modules (`emoteInput`, `emoteStat`, `checkEmotes`, `detTimeMem`). Step II is performed by `emote_step_II.py` (this repository) instead of the original `emoteStepII.pl`

The paths to these programs are hardcoded in the script (see above). `requests` is used to retrieve pathway information from the KEGG database for pathway enrichment analysis, so internet access is required when running `analysis.py`.

## License

Released under the MIT License.

```
MIT License

Copyright (c) 2026 Dr. Knut Finstermeier, Max Planck Unit for the Science of Pathogens

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
