# plaseek

構造ベースでプラスミド内のタンパク質を検索するツール

## 必要なもの

### バイナリ

- [Python](https://www.python.org/) 3.10 or later
- [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [GNU parallel](https://www.gnu.org/software/parallel/)
- [Foldseek](https://github.com/steineggerlab/foldseek)

これらのバイナリは[Homebrew](https://brew.sh/)を使ってインストールできます。Homebrewをインストールした後、ターミナルで以下のコマンドを実行してください。

```bash
brew install python@3.11 blast parallel brewsci/bio/foldseek
```

### データベース

`makeblastdb`を使ってターゲット配列データベースを準備しておきます。

```bash
makeblastdb -in Pseudomonas_plasmids_847.fasta -dbtype nucl -out P847DB -parse_seqids
```

### インストール

```bash
# install from GitHub
python3.11 -m pip install git+https://github.com/YoshitakaMo/plaseek.git
# upgrade
python3.11 -m pip uninstall plaseek -y && python3.11 -m pip install --upgrade git+https://github.com/YoshitakaMo/plaseek.git
```

## plaseekの使い方

インプットファイルとしてPDB(`.pdb`)ファイルまたは[Foldseek Search Server](https://search.foldseek.com/search)から得られたm8ファイル(`.m8`)を指定します。インプットファイルは`-i`オプションで指定します。

- 引数`-i`には入力ファイルを指定します。**設定は必須です**。`-i`には`.pdb`または`.m8`ファイルのいずれか一方を指定できます。PDBファイルを指定した場合、かつローカルに存在するFoldseekのバイナリ（`--foldseek_binary_path`）とデータベース（`--foldseek_db_path`または`-f`）へのPATHが未設定だった場合、自動的にインターネット上の[Foldseek Search Server](https://search.foldseek.com/search)を使って類縁構造の結果ファイルを取得した後、BLASTを使ってターゲット配列データベースに対して検索を行います。m8ファイルを指定した場合はFoldseekを使った類縁構造検索は省略されます。
- 引数`-t`または`--target-sequence-db-path`には検索対象とするターゲット配列データベースへのPATHを設定します。**設定は必須です**。設定方法は上記の[データベース](#データベース)の項を参照してください。
- 引数`-o`に結果ファイルの出力先ディレクトリを指定します。未入力の場合は`-i`で入力したファイル名から拡張子を除いた名前のディレクトリに出力されます。
-

```bash
# show help message
plaseek -h

# run plaseek for (predicted) pdb file
# Use Foldseek web server to search for homologous structures if both `--foldseek_binary_path` and `--foldseek_db_path (-f)` are not specified.
plaseek -i AF-P07676-F1-model_v4.pdb \
        -t /path/to/db/P847DB

# run plaseek for (predicted) pdb file with locally-installed Foldseek
# The result files will be saved in `outputresults` specified by the `-o` arg.
plaseek -i AF-P07676-F1-model_v4.pdb \
        -t /path/to/db/P847DB \
        --foldseek_binary_path /home/moriwaki/apps/bin/foldseek \
        --foldseek_db_path /scr/foldseekdb
        -o outputresults

# run plaseek for m8 file obtained from Foldseek web server
# The result files will be saved in `result2` directory.
plaseek -i alis_afdb50.m8 \
        -t /path/to/db/P847DB \
        -o result2
```

`--foldseek_evalue_threshold`と`--tblastn_pident_threshold`の値を調節することで検索結果の精度を変更することができます。`--foldseek_evalue_threshold`の値を大きくすると、より広い範囲の類縁構造を検索することができますが、検索結果の精度は下がります。`--tblastn_pident_threshold`の値を小さくすると、target sequence DBの中からその一致度以下のアミノ酸配列を持つ配列を取得することができますが、精度が低下します。

```bash
# Foldseek e-value threshold: 1e-10 (default: 1e-20), tblastn pident threshold: 90.0% (default: 98.0%)
plaseek -i alis_afdb50.m8 \
        -t /path/to/db/P847DB \
        -o result3 \
        --foldseek_evalue_threshold 1e-10 \
        --tblastn_pident_threshold 90.0
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
