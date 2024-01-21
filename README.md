# repseek

Search for plamid replication proteins based on protein structures

## prerequisites

### binaries

- [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [GNU parallel](https://www.gnu.org/software/parallel/)
- [Foldseek](https://github.com/steineggerlab/foldseek)

These binaries can be installed using [Homebrew](https://brew.sh/). After installing Homebrew, run the following commands in a terminal:

```bash
brew install blast parallel
brew install brewsci/bio/foldseek
```

### databases

It is recommended to use `makeblastdb` to prepare target sequence databases:

```bash
makeblastdb -in Pseudomonas_plasmids_847.fasta -dbtype nucl -out P847DB -parse_seqids
```

## Usage

```bash
# show help message
repseek -h
# run repseek with default parameters
repseek -i AF-P07676-F1-model_v4.pdb \
        --foldseek-db-path /scr/foldseek \
        --target-sequence-db-path ../../db/P847DB \
        -o results.txt
```