# plaseek

Search for plamid replication proteins based on protein structures

## prerequisites

### binaries

- [Python](https://www.python.org/) 3.9 or later
- [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [GNU parallel](https://www.gnu.org/software/parallel/)
- [Foldseek](https://github.com/steineggerlab/foldseek)

These binaries can be installed using [Homebrew](https://brew.sh/). After installing Homebrew, run the following commands in a terminal:

```bash
brew install blast parallel brewsci/bio/foldseek
```

### databases

Use `makeblastdb` to prepare target sequence databases:

```bash
makeblastdb -in Pseudomonas_plasmids_847.fasta -dbtype nucl -out P847DB -parse_seqids
```

### installation

```bash
# install from GitHub
python3.12 -m pip install git+https://github.com/YoshitakaMo/plaseek.git
# upgrade
python3.12 -m pip install --upgrade git+https://github.com/YoshitakaMo/plaseek.git
```

## Usage

```bash
# show help message
plaseek -h

# run plaseek for (predicted) pdb file
plaseek -i AF-P07676-F1-model_v4.pdb \
        --foldseek-db-path /scr/foldseek \
        --target-sequence-db-path ../../db/P847DB \
        -o results.txt

# run plaseek for m8 file obtained from Foldseek web server
plaseek -i AF-P07676-F1-model_v4.m8 \
        --foldseek-db-path /scr/foldseek \
        --target-sequence-db-path ../../db/P847DB \
        -o results.txt
```