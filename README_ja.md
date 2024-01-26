# plaseek

構造ベースでプラスミド内の複製タンパク質を検索するツール

## 必要なもの

### バイナリ

- [Python](https://www.python.org/) 3.9 or later
- [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [GNU parallel](https://www.gnu.org/software/parallel/)
- [Foldseek](https://github.com/steineggerlab/foldseek)

これらのバイナリは[Homebrew](https://brew.sh/)を使ってインストールできます。Homebrewをインストールした後、ターミナルで以下のコマンドを実行してください。

```bash
brew install python@3.12 blast parallel brewsci/bio/foldseek
```

### データベース

`makeblastdb`を使ってターゲット配列データベースを準備しておきます。

```bash
makeblastdb -in Pseudomonas_plasmids_847.fasta -dbtype nucl -out P847DB -parse_seqids
```

### インストール

```bash
# install from GitHub
python3.12 -m pip install git+https://github.com/YoshitakaMo/plaseek.git
# upgrade
python3.12 -m pip install --upgrade git+https://github.com/YoshitakaMo/plaseek.git
```

## 使い方

インプットファイルとしてPDB(`.pdb`)ファイルまたはFoldseekのウェブサーバーから得られたm8ファイル(`.m8`)を指定します。インプットファイルは`-i`オプションで指定します。

PDBファイルを指定した場合、Foldseekを使って類縁構造を検索した後、BLASTを使ってターゲット配列データベースに対して検索を行います。m8ファイルを指定した場合はFoldseekを使った類縁構造検索は省略されます。

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
