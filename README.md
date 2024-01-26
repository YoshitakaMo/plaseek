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

## plaseekの使い方

インプットファイルとしてPDB(`.pdb`)ファイルまたは[Foldseek Search Server](https://search.foldseek.com/search)から得られたm8ファイル(`.m8`)を指定します。インプットファイルは`-i`オプションで指定します。

PDBファイルを指定した場合、Foldseekを使って類縁構造を検索した後、BLASTを使ってターゲット配列データベースに対して検索を行います。m8ファイルを指定した場合はFoldseekを使った類縁構造検索は省略されます。

```bash
# show help message
plaseek -h

# run plaseek for (predicted) pdb file
plaseek -i AF-P07676-F1-model_v4.pdb \
        --foldseek-db-path /scr/foldseek \
        --target-sequence-db-path /path/to/db/P847DB \
        -o results.txt

# run plaseek for m8 file obtained from Foldseek web server
plaseek -i AF-P07676-F1-model_v4.m8 \
        --target-sequence-db-path /path/to/db/P847DB \
        -o results.txt
```

## Foldseek Search Serverで類縁構造のリストを作成する

Foldseek Search Serverで類縁構造のリストを作成するには、以下の手順を実行します。

### Foldseek Search Serverで検索を行う

検索したい配列のUniProt accession IDがわかっている場合は、LOAD ACESSIONボタンを押してそのAccession IDを入力すると、AlphaFold2予測構造が自動的にロードされます。または、UPLOAD PDBボタンから手持ちのPDBファイルをアップロードすることができます。

![upload pdb](https://i.imgur.com/nGGYL6t.png)

### Databases & search settingsのパラメータ設定

DatabasesはAlphaFold/UniProt50(v4)にチェックを入れておきます。後はPDB100も含めて通常は外しておいて構いませんが、MGnifyなどは付け加えてみても良いかもしれません。Modeは3Di/AAの方が早く検索結果が帰ってきます。
Taxnomic filterを使うと、検索結果を特定の分類群に限定することができます。

![databases](https://i.imgur.com/3jqb4ze.png)

### Searchを実行する

SEARCHボタンを押して1〜2分待つと結果が表示されます。

![search](https://i.imgur.com/xneRSS9.png)

左上にメニューボタンがあり、そこから検索結果をダウンロードすることができます。

![search results](https://i.imgur.com/xGxZML5.png)

Hit tablesからM8ファイル（が圧縮されたtar.gzファイル）をダウンロードします。

![hit tables](https://i.imgur.com/4s6uMkX.png)

ファイルを解凍すると、plaseekに用いることのできるM8ファイルの構成が得られます。複数のデータベースを指定した場合はその数だけのm8ファイルが得られます。

![m8files](https://i.imgur.com/hgVDmu6.png)

