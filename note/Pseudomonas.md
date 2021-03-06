# *Pseudomonas* HGT

This markdown is to record my study of teacher's [Pseudomonas.md](https://github.com/wang-q/withncbi/blob/master/pop/Pseudomonas.md).

Before getting start, there is a little background information for better understanding the following processes.

## Background info

- **Taxonomy**

The basic taxonomy ranks are here:

| Domain  | Kingdom  | Phylum | Class   | Order  | Family   | Genus  | Species |
|---------|----------|--------|---------|--------|----------|--------|---------|
| Domains | Kingdoms | Phyla  | Classes | Orders | Families | Genera | Species |
| 域       | 界        | 门      | 纲       | 目      | 科        | 属      | 种       |

**Knowledge Point**:

[Back-formation](https://en.wiktionary.org/wiki/Appendix:Glossary#back-formation): A term formed by removing an apparent or real prefix or suffix from an older term; just shortens a word without changing.

Taxon is back-formation from taxonomy. And its plural format is taxa.

- **Prokaryotic RefSeq Genomes**

[Prokaryotic RefSeq Genomes](https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/) included archaeal and bacterial genome assemblies which meet sequence and annotation quality criteria. There are different level for you to understand.

1. RefSeq: as mentioned above, included genome assemblies that meet quality criteria.
2. Reference genomes: genome assemblies that are annotated and updated by the assembly submitters and chosen by the RefSeq curatorial staff based on their quality and importance to the community as anchors for the analysis of other genomes in their taxonomic group. Reference genomes are annotated with YP_ or NP_protein accessions.
3. Representative genomes: for species without a reference genome, one assembly per defined species is selected as representative.

- **Markdown goal**:

The genus Pseudomonas includes the conditionally pathogenic bacteria Pseudomonas aeruginosa, plant pathogens, plant beneficial bacteria, and soil bacteria. Microorganisms of Pseudomonas are extremely rich in metabolic diversity, and it is thought that this diversity also allows them to survive in a very wide range of ecological niches.

The metabolism of carbohydrates is a fundamental biochemical process that ensures a continuous supply of energy to living cells.

The biodegradation of RNA is also an important part of metabolism, which is accomplished by the degradosomes. However, the diversity of degradosomes in different environments has not been fully investigated.

According to a recent [paper](https://journals.asm.org/doi/10.1128/mSystems.00543-20), there are some order-level changes in Gammaproteobacteria. We include both old and new orders.

Old ones: Cellvibrionales, Oceanospirillales, Pseudomonadales, and Alteromonadales
New ones: Moraxellales, Kangiellales, and Pseudomonadales

## Software

`PPanGGOLiN` and `InterProScan` installation steps can be found [here](https://github.com/wang-q/dotfiles/blob/master/others.sh)

- PPanGGOLiN

Need python_3.9

```bash
if grep -q -i PYTHON_39_PATH $HOME/.bashrc; then
    echo "==> .bashrc already contains PYTHON_39_PATH"
else
    echo "==> Update .bashrc"

    echo '# PYTHON_39_PATH' >> $HOME/.bashrc
    PYTHON_39_PATH="export PATH=\"/home/linuxbrew/.linuxbrew/Cellar/python@3.9/3.9.13/bin:\$PATH\""
    echo ${PYTHON_39_PATH} >> $HOME/.bashrc
fi
```

- InterProScan

Go check [others.sh](https://github.com/wang-q/dotfiles/blob/master/others.sh)

- `nwr`

```bash
brew install wang-q/tap/nwr # 0.5.5 or above
brew install sqlite         # 3.34 or above

nwr download
nwr txdb

nwr ardb
nwr ardb --genbank
```

- Other packages

If you already install those packages by `brew`, then nothing will happen.

So better check them whether you have installed them or not.

```bash
brew install hmmer
brew install brewsci/bio/muscle
brew install brewsci/bio/fasttree
brew install brewsci/bio/newick-utils
brew install brewsci/bio/trimal

brew install datamash
brew install miller
brew install wang-q/tap/tsv-utils

brew install librsvg
brew install jq
brew install pup

# dN/dS
brew install brewsci/bio/clustal-w
brew install brewsci/bio/paml

cpanm Bio::Tools::Run::Alignment::Clustalw
cpanm https://github.com/wang-q/Bio-Tools-Phylo-PAML.git
```

- Tree viewing tools on windows

```powershell
winget install -s winget -e --id launch4j.launch4j
winget install -s winget -e --id Oracle.JavaRuntimeEnvironment
```

```bash
cd /mnt/d/download
wget https://github.com/rambaut/figtree/releases/download/v1.4.4/FigTree_v1.4.4.tgz

cd /mnt/d/Biosoftware
tar -xzvf ../download/FigTree_v1.4.4.tgz
```

If you do not have JavaRuntimeEnvironment, the `FigTree` will not run properly.

The `.jar` file in `./figtree_v1.4.4/lib` is the right one.

## Strain info

Genera:

- Pseudomonas
- Acinetobacter

### List all ranks

```bash
mkdir -p /mnt/e/data/Pseudomonas
cd /mnt/e/data/Pseudomonas

# count Pseudomonas by nwr
nwr member Pseudomonas |
    grep -v " sp." |
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat

# count Acinetobacter by nwr
nwr member Acinetobacter |
    grep -v " sp." |
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat

# count species and species subgroup in both Genus
nwr member Pseudomonas Acinetobacter -r "species group" -r "species subgroup" |
    tsv-select -f 1-3 |
    keep-header -- tsv-sort -k3,3 -k2,2 |
    sed 's/Pseudomonas /P. /g' |
    sed 's/Acinetobacter /A. /g' |
    mlr --itsv --omd cat
```

Pseudomonas:

| rank             | count |
|------------------|-------|
| genus            | 1     |
| species          | 415   |
| strain           | 747   |
| subspecies       | 12    |
| no rank          | 120   |
| species group    | 6     |
| species subgroup | 5     |
| isolate          | 1     |

Acinetobacter:

| rank             | count |
|------------------|-------|
| genus            | 1     |
| species group    | 2     |
| species subgroup | 3     |
| species          | 114   |
| strain           | 1110  |
| no rank          | 2     |
| subspecies       | 1     |
| isolate          | 2     |

Species group/subgroup:

| #tax_id | sci_name                                 | rank             |
|---------|------------------------------------------|------------------|
| 2839056 | A. Taxon 24                              | species group    |
| 909768  | A. calcoaceticus/baumannii complex       | species group    |
| 136841  | P. aeruginosa group                      | species group    |
| 136842  | P. chlororaphis group                    | species group    |
| 136843  | P. fluorescens group                     | species group    |
| 136845  | P. putida group                          | species group    |
| 136846  | P. stutzeri group                        | species group    |
| 136849  | P. syringae group                        | species group    |
| 2839060 | A. Taxon 24C                             | species subgroup |
| 2839057 | A. Taxon 24D                             | species subgroup |
| 2839061 | A. Taxon 24E                             | species subgroup |
| 627141  | P. nitroreducens/multiresinivorans group | species subgroup |
| 1232139 | P. oleovorans/pseudoalcaligenes group    | species subgroup |
| 578833  | P. stutzeri subgroup                     | species subgroup |
| 251695  | P. syringae group genomosp. 1            | species subgroup |
| 251698  | P. syringae group genomosp. 2            | species subgroup |

### Species with assemblies

Also check the order Pseudomonadales.

- Old ones: Cellvibrionales, Oceanospirillales, Pseudomonadales, and Alteromonadales
- New ones: Moraxellales, Kangiellales, and Pseudomonadales

Here means that orders have been changed. See [above](#pseudomonas-hgt) for the paper offered.
  
```bash
cd /mnt/e/data/Pseudomonas

SPECIES=$(
    nwr member -r species \
        Cellvibrionales Oceanospirillales Alteromonadales \
        Moraxellales Kangiellales Pseudomonadales |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        grep -v " sp." |
        grep -v -E "\bbacterium\b" |
        grep -v -E "\bsymbiont\b" |
        sed '1d' |
        cut -f 1 |
        sort |
        uniq
)

for S in $SPECIES; do
    RS=$(
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND species_id = $S
            " |
            sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
    )

    CHR=$(
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND species_id = $S
                AND assembly_level IN ('Complete Genome', 'Chromosome')
            " |
            sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
    )

    if [[ ${RS} -gt 0 ]]; then
        echo -e "$S\t$RS\t$CHR"
    fi
done |
    nwr append stdin |
    tsv-select -f 1,4,2-3 |
    tsv-sort -k3,3nr -k4,4nr -k2,2 |
    (echo -e 'species_id\tspecies\tRS\tCHR' && cat) \
    > species.count.tsv

# strains have more than 5 assemblies in two genera
cat species.count.tsv |
    tsv-filter -H --ge CHR:5 |
    tsv-filter -H --invert --str-in-fld species:Pseudomonas --lt RS:30 |
    tsv-filter -H --invert --str-in-fld species:Acinetobacter --lt RS:30 |
    sed 's/Pseudomonas /P. /g' |
    sed 's/Acinetobacter /A. /g' |
    mlr --itsv --omd cat

# grep:
# -i/--ignore-case: force grep to ignore word case (boo, Boo, BOO, ...)
# -E/--extended-regexp: Treats pattern as an extended regular expression (ERE)
# \b allows you to perform a "whole words only" search using a regular expression in the form of \bword\b

# nwr:
# nwr member:
# output format: #tax_id sci_name        rank    division
# nwr append:
# append fields of higher ranks to a TSV file
# * If `--rank` is empty, the scientific name will be appended.

# tsv-sort:
# tsv-sort runs 'sort' using TAB as the field delimiter, all arguments are forwarded to 'sort'
# -r/--reverse: reverse the result of comparisons
# -n/--numeric-sort: compare according to string numerical value
# -k/--key=KEYDEF: sort via a key; KEYDEF - pos1[,pos2]
# KEYDEF: Specify a sort field that consists of the part of the line between pos1 and pos2 
#         (or the end of the line, if pos2 is omitted), inclusive.
#         fields being separated by runs of blank characters

# (echo -e "" && cat) - means echo -e and cat contents from stdin ("|")
# so the step will combined echo -e and stdin
```

> **Explanation**:
>
> The first step was acquiring all tax_ids from col1 into `$SPECIES` as an array.
>
> ```txt
> sqlite3: DB-API 2.0 interface for SQLite databases
> FILENAME is the name of an SQLite database. A new database is created if the
> file does not previously exist.
> 
> Usage: sqlite3 [OPTIONS] FILENAME [SQL]
> 
> -tabs: set output mode to 'tabs'
> ```
>
> `[SQL]` means the SQLite language
>
> ```txt
> # select elements according to conditions combined
> SELECT column1, column2....columnN
> FROM   table_name
> WHERE  CONDITION-1 {AND|OR} CONDITION-2;
>
> # count all elements selected
> SELECT COUNT(column_name)
> FROM   table_name
> WHERE  CONDITION;
> ```
>
> So the `$RS` were numbers of all raws counted from a species_id
>
> And the `$CHR` were numbers of all raws counted from a species_id which assembly level were "Complete Genome" and "Chromosome".

**Notice**:

Basic (BRE) and extended (ERE) regular expression are two variations on the syntax of the specified pattern. BRE is the default in `sed` and `grep`. Use the POSIX-specified `-E` option (`-r, --regexp-extended`) to enable Extended Regular Expression (ERE) syntax. (The contents were from [BRE-vs-ERE](https://www.gnu.org/software/sed/manual/html_node/BRE-vs-ERE.html))

| species_id | species                         | RS   | CHR |
|------------|---------------------------------|------|-----|
| 287        | P. aeruginosa                   | 6607 | 533 |
| 470        | A. baumannii                    | 6340 | 389 |
| 33069      | P. viridiflava                  | 1538 | 7   |
| 317        | P. syringae                     | 598  | 45  |
| 48296      | A. pittii                       | 330  | 32  |
| 294        | P. fluorescens                  | 259  | 39  |
| 480        | Moraxella catarrhalis           | 209  | 16  |
| 303        | P. putida                       | 194  | 50  |
| 106654     | A. nosocomialis                 | 165  | 11  |
| 316        | P. stutzeri                     | 135  | 31  |
| 29438      | P. savastanoi                   | 116  | 5   |
| 38313      | Shewanella algae                | 109  | 22  |
| 587753     | P. chlororaphis                 | 99   | 60  |
| 756892     | A. indicus                      | 93   | 21  |
| 47877      | P. amygdali                     | 89   | 10  |
| 380021     | P. protegens                    | 73   | 23  |
| 40215      | A. junii                        | 66   | 10  |
| 1530123    | A. seifertii                    | 60   | 25  |
| 29430      | A. haemolyticus                 | 55   | 14  |
| 76759      | P. monteilii                    | 47   | 9   |
| 40214      | A. johnsonii                    | 43   | 19  |
| 40216      | A. radioresistens               | 42   | 5   |
| 296        | P. fragi                        | 39   | 6   |
| 28090      | A. lwoffii                      | 33   | 11  |
| 43657      | Pseudoalteromonas luteoviolacea | 25   | 5   |
| 28108      | Alteromonas macleodii           | 24   | 9   |
| 24         | Shewanella putrefaciens         | 23   | 12  |
| 34062      | Moraxella osloensis             | 23   | 10  |
| 43662      | Pseudoalteromonas piscicida     | 19   | 6   |
| 314275     | Alteromonas mediterranea        | 16   | 16  |
| 62322      | Shewanella baltica              | 14   | 11  |
| 386891     | Moraxella bovoculi              | 9    | 7   |
| 1697053    | Thiopseudomonas alkaliphila     | 7    | 7   |

### Outgroups

Use these model organisms as outgroups.

```bash
cd /mnt/e/data/Pseudomonas

# other bacteria genera
GENUS=$(
    nwr member Bacteria -r genus |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        sed '1d' |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//' # remove the last ','
)

# extract reference genome from all genera
echo "
.headers ON

    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus_id IN ($GENUS)
        AND refseq_category IN ('reference genome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > reference.tsv

cat reference.tsv | wc -l
#  16 - only 16 reference genome from all genera

cat reference.tsv |
    sed '1s/^/#/' |
    nwr append stdin -r phylum -r class |
    tsv-select -H -f 1,2,phylum,class |
    parallel --col-sep "\t" -j 1 '
        if [[ "{3}" != "Proteobacteria" ]]; then
            printf "%s\t%s\t%s\n" {1} {2} {3}
        else
            printf "%s\t%s\t%s\n" {1} {2} {3}/{4}
        fi
    ' |
    mlr --itsv --omd cat

# sed: add # to the start of first line
# nwr append: append phylum and class to the tsv
```

Reference genome:

| #tax_id | organism_name                                                    | phylum                               |
|---------|------------------------------------------------------------------|--------------------------------------|
| 565050  | Caulobacter vibrioides NA1000                                    | Proteobacteria/Alphaproteobacteria   |
| 192222  | Campylobacter jejuni subsp. jejuni NCTC 11168 = ATCC 700819      | Proteobacteria/Epsilonproteobacteria |
| 208964  | Pseudomonas aeruginosa PAO1                                      | Proteobacteria/Gammaproteobacteria   |
| 871585  | Acinetobacter pittii PHEA-2                                      | Proteobacteria/Gammaproteobacteria   |
| 511145  | Escherichia coli str. K-12 substr. MG1655                        | Proteobacteria/Gammaproteobacteria   |
| 386585  | Escherichia coli O157:H7 str. Sakai                              | Proteobacteria/Gammaproteobacteria   |
| 1125630 | Klebsiella pneumoniae subsp. pneumoniae HS11286                  | Proteobacteria/Gammaproteobacteria   |
| 99287   | Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 | Proteobacteria/Gammaproteobacteria   |
| 198214  | Shigella flexneri 2a str. 301                                    | Proteobacteria/Gammaproteobacteria   |
| 227377  | Coxiella burnetii RSA 493                                        | Proteobacteria/Gammaproteobacteria   |
| 272561  | Chlamydia trachomatis D/UW-3/CX                                  | Chlamydiae                           |
| 93061   | Staphylococcus aureus subsp. aureus NCTC 8325                    | Firmicutes                           |
| 224308  | Bacillus subtilis subsp. subtilis str. 168                       | Firmicutes                           |
| 169963  | Listeria monocytogenes EGD-e                                     | Firmicutes                           |
| 83332   | Mycobacterium tuberculosis H37Rv                                 | Actinobacteria                       |

Actinobacteria, Firmicutes are phyla of Gram-positive bacteria.

Proteobacteria(synonym Pseudomonadota) is a major phylum of Gram-negative bacteria.

Chalamydiae are obligate intracellular parasites.



## Download all assemblies

- Get all assemblies info

```bash
cd /mnt/e/data/Pseudomonas

# Phylum
# Pseudomonas aeruginosa PAO1 is in the reference list
# All Refseqs from different phyla were contained as outgroup
cat reference.tsv |
    tsv-select -H -f organism_name,species,genus,ftp_path,assembly_level \
    > raw.tsv
cat raw.tsv | wc -l
#  16

# Genus (contained all Species)
# Species with 2 or more genomes were retained from our target genera (changed order)
SPECIES=$(
    cat species.count.tsv |
        tsv-filter -H --ge CHR:2 |
        sed '1d' |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)
# because of the CHR meaning numbers of complete genome/chromosome inside a species
# so using CHR >= 2 to keep those species

# Order
# Pseudomonadales
echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    tsv-filter --invert --str-eq 2:"Pseudomonas aeruginosa" --str-eq 5:"Chromosome" |
    tsv-filter --invert --str-eq 2:"Acinetobacter baumannii" --str-eq 5:"Chromosome" \
    >> raw.tsv
# Other strains belong to Pseudomonadales order but not our target species

cat raw.tsv | wc -l
#  1666

# The SQLite || operator allows you to concatenate 2 or more strings together.
# This step finally got all strains' assemblies from species.count.tsv

# Class
# Also includes representative strains of Gammaproteobacteria (Class).
# Representative strains were unequal to RefSeq
# families not in our orders
FAMILY=$(
    nwr member Gammaproteobacteria -r family |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        sed '1d' |
        cut -f 1 |
        nwr append stdin -r order |
        tsv-filter --str-ne 2:"Cellvibrionales" --str-ne 2:"Oceanospirillales" --str-ne 2:"Alteromonadales" |
        tsv-filter --str-ne 2:"Moraxellales" --str-ne 2:"Kangiellales" --str-ne 2:"Pseudomonadales" |
        tsv-filter --str-ne 2:"NA" |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND family_id IN ($FAMILY)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND refseq_category IN ('representative genome')
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    grep -v -i "symbiont " |
    tsv-filter --str-not-in-fld 1:"[" \
    >> raw.tsv

cat raw.tsv | wc -l
#  2095

# Create abbr.
cat raw.tsv |
    grep -v '^#' |
    tsv-uniq |
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "1,2,3" -s '\t' -m 3 --shortsub |
    (echo -e '#name\tftp_path\torganism\tassembly_level' && cat ) |
    perl -nl -a -F"," -e '
        BEGIN{my %seen};
        /^#/ and print and next;
        /^organism_name/i and next;
        $seen{$F[3]}++; # ftp_path
        $seen{$F[3]} > 1 and next;
        $seen{$F[5]}++; # abbr_name
        $seen{$F[5]} > 1 and next;
        printf qq{%s\t%s\t%s\t%s\n}, $F[5], $F[3], $F[1], $F[4];
        ' |
    keep-header -- sort -k3,3 -k1,1 \
    > Pseudomonas.assembly.tsv
# %seen in perl means only print the first line that was seen in the input
# A double dash (--) delimits the command in keep-header, similar to how the pipe operator (|) delimits commands. Examples:
# keep-header file1.txt -- sort

# find potential duplicate strains or assemblies
cat Pseudomonas.assembly.tsv |
    tsv-uniq -f 1 --repeated
# --repeated: output only lines that are repeated (based on the key)

# Edit .tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Pseudomonas.assembly.tsv
# cp Pseudomonas.assembly.tsv ~/Scripts/withncbi/pop

# Comment out unneeded strains

# Cleaning
rm raw*.*sv
```

> **Explanation**:
>
> [SQLite: || Operator](https://www.techonthenet.com/sqlite/functions/concatenate.php)
>
> `||` operator allows you to concatenate 2 or more strings together.
>
> ```txt
> SELECT 'Jane' || ' ' || 'Smith';
> Result: 'Jane Smith'
> ```
>
> The `||` operator will concatenate string values that are enclosed in single quotes.
>
> `LIKE '% sp.%'` clause can find which tracks composed by `sp.`
>
> `%` here represents other characters:
>
> `'% sp.'` will match those tracks ended with `sp.`, such as `abcsp.`
>
> `'sp.%'` will match those tracks started with `sp.`, such as `sp.abc`
>
> `'% sp.%'` will match those tracks contained `sp.`
>
> ```bash
> perl ~/Scripts/withncbi/taxon/abbr_name.pl
> Usage:
>       cat <file> | perl abbr_name.pl [options]
>         Options:
>           --help              brief help message
>           --column    -c  STR Columns of strain, species, genus, default is 1,2,3.
>                               If there's no strain, use 1,1,2.
>                               Don't need the strain part, use 2,2,3
>                               When there's only strain, use 1,1,1
>           --separator -s  STR separator of the line, default is "\s+"
>           --min INT           mininal length for abbreviation of species
>           --tight             no underscore between Genus and species
>           --shortsub          clean subspecies parts
> ```

- Download them using scripts

```bash
cd /mnt/e/data/Pseudomonas

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f Pseudomonas.assembly.tsv \
    -o ASSEMBLY

# Remove dirs not in the list
find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
    tr "/" "\t" |
    cut -f 2 |
    tsv-join --exclude -k 1 -f ASSEMBLY/rsync.tsv -d 1 |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        echo Remove {}
        rm -fr ASSEMBLY/{}
    '

# Run
proxychains4 bash ASSEMBLY/Pseudomonas.assembly.rsync.sh

bash ASSEMBLY/Pseudomonas.assembly.collect.sh

# md5
cat ASSEMBLY/rsync.tsv |
    tsv-select -f 1 |
    parallel -j 4 --keep-order '
        echo "==> {}"
        cd ASSEMBLY/{}
        md5sum --check md5checksums.txt
    ' |
    grep -v ": OK"
```

> **Explanation**:
>
> ```bash
> perl ~/Scripts/withncbi/taxon/assembly_prep.pl --help
> Usage:
>       perl assembly_prep.pl [options]
>         Options:
>           --help, -?              brief help message
>
>           --file, -f      STR     tab seperated file containing wgs prefix and name
>           --outdir, -o    STR     output dir
>           --csvonly               only generate the csv file
>
>       perl assembly_prep.pl -f trichoderma.assembly.tsv -o ASSEMBLY
>
>       #name   ftp_path    organism    assembly_level
>
>       Three files will be generated:
>
>           trichoderma.assembly.rsync.sh
>           rsync.tsv
>           trichoderma.assembly.collect.sh
>
>       The latter one will create:
>
>           trichoderma.assembly.collect.csv
> ```

## BioSample

- Extract all samples' info into a tsv

ENA's BioSample missed many strains, so NCBI's was used.

```bash
cd /mnt/e/data/Pseudomonas

mkdir -p biosample

ulimit -n `ulimit -Hn`

cat ASSEMBLY/Pseudomonas.assembly.collect.csv | 
    tsv-select -H -d, -f BioSample | 
    head
#BioSample
#SAMN05526184
#SAMN05853496
#SAMN06240324
#SAMN05828143
#SAMN02603494
#SAMN02603051
#SAMN02603889
#SAMN02603140
#SAMN09302593

cat ASSEMBLY/Pseudomonas.assembly.collect.csv |
    tsv-select -H -d, -f BioSample |
    grep "^SAM" |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        if [ ! -s biosample/{}.txt ]; then
            >&2 echo {}
            curl -fsSL "https://www.ncbi.nlm.nih.gov/biosample/?term={}&report=full&format=text" -o biosample/{}.txt
#            curl -fsSL "https://www.ebi.ac.uk/biosamples/samples/{}" -o biosample/{}.json
        fi
    '
# curl:
# -f/--fail: fail fast with no output on HTTP errors
# -s/--silent: silent mode
# -S/--show-error: show error even when -s is used
# -L/--location: follow redirects
# This step will actually give you all the Biosample into a txt according to their SAMN number

find biosample -name "SAM*.txt" | wc -l
#  1957

# get all BioSample txt attributes
find biosample -name "SAM*.txt" |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        cat {} |
            perl -nl -e '\''
                print $1 if m{\s+\/([\w_ ]+)=};
            '\''
    ' |
    tsv-uniq --at-least 50 | # ignore rare attributes
    grep -v "^INSDC" |
    grep -v "^ENA" \
    > attributes.lst

cat attributes.lst |
    (echo -e "BioSample" && cat) |
    tr '\n' '\t' |
    sed 's/\t$/\n/' \
    > Pseudomonas.biosample.tsv
# tr: transform line to header splited by tab
# sed: change the tab at the line end to return
# make attributes.lst a line as header of Pseudomonas.biosample.tsv

find biosample -name "SAM*.txt" |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        >&2 echo {/.}
        cat {} |
            perl -nl -MPath::Tiny -e '\''
                BEGIN {
                    our @keys = grep {/\S/} path(q{attributes.lst})->lines({chomp => 1});
                    our %stat = ();
                }

                m(\s+\/([\w_ ]+)=\"(.+)\") or next;
                my $k = $1;
                my $v = $2;
                if ( $v =~ m(\bNA|missing|Not applicable|not collected|not available|not provided|N\/A|not known|unknown\b)i ) {
                    $stat{$k} = q();
                } else {
                    $stat{$k} = $v;
                }

                END {
                    my @c;
                    for my $key ( @keys ) {
                        if (exists $stat{$key}) {
                            push @c, $stat{$key};
                        }
                        else {
                            push @c, q();
                        }
                    }
                    print join(qq{\t}, q{{/.}}, @c);
                }
            '\''
    ' \
    >> Pseudomonas.biosample.tsv
# This script actually extract each sample attributes into a whole tsv
# @keys are the target attributes in the headline
```

> ```perl
> # reading files
> @lines = $file->lines;
> @contents = path("/tmp/foo.txt")->lines( { chomp => 1, count => 4 } );
> ```

## Count and group strains

- Check N50 of assemblies

- Some strains were anomalously labeled and identified by the `mash` tree.
  - Pseudom_flu_GCF_900636635_1
  - Pseudom_chl_GCF_001023535_1
  - Pseudom_syr_GCF_004006335_1
  - Pseudom_puti_GCF_003228315_1 and Pseudom_puti_GCF_020172705_1

Those strains have been removed after [Raw phylogenetic tree by MinHash](#raw-phylogenetic-tree-by-minhash) completed. So actually this step was repeated several times until all abnormal strains removed properly. Those strains were not removed at the outset, but underwent repeated following validation from here to [Raw phylogenetic tree by MinHash](#raw-phylogenetic-tree-by-minhash).

```bash
cd /mnt/e/data/Pseudomonas

for dir in $(find ASSEMBLY -maxdepth 1 -mindepth 1 -type d | sort); do
    1>&2 echo "==> ${dir}"
    name=$(basename ${dir})

    find ${dir} -type f -name "*_genomic.fna.gz" |
        grep -v "_from_" | # exclude CDS and rna
        xargs cat |
        faops n50 -C -S stdin |
        (echo -e "name\t${name}" && cat) |
        datamash transpose
done |
    tsv-uniq | # remove all header, keep the first line
    tee ASSEMBLY/n50.tsv
# n50.tsv after transpose will give out every name N50 info
# xargs cat: because the faops accepted the gz file, so cat fna.gz and xargs to faops
# file header: name\tN50\tS\tC
# -S: compute sum of size of all entries
# -C: count entries

cat ASSEMBLY/n50.tsv |
    tsv-filter \
        -H --or \
        --le 4:100 \
        --ge 2:100000 |
    tsv-filter -H --ge 3:1000000 |
    tr "\t" "," \
    > ASSEMBLY/n50.pass.csv
# evaluate tests as an OR rather than an AND clause (default)

wc -l ASSEMBLY/n50*
#  1958 ASSEMBLY/n50.pass.csv
#  1959 ASSEMBLY/n50.tsv
#  3917 total

tsv-join \
    ASSEMBLY/Pseudomonas.assembly.collect.csv \
    --delimiter "," -H --key-fields 1 \
    --filter-file ASSEMBLY/n50.pass.csv |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_900636635 |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_001023535 |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_004006335 |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_003228315 |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_020172705 \
    > ASSEMBLY/Pseudomonas.assembly.pass.csv

wc -l ASSEMBLY/Pseudomonas.assembly*csv
#  1959 ASSEMBLY/Pseudomonas.assembly.collect.csv
#  1953 ASSEMBLY/Pseudomonas.assembly.pass.csv

# number of strains
cat ASSEMBLY/Pseudomonas.assembly.pass.csv | sed '1d' | wc -l
#  1952
```

- Order

```bash
cd /mnt/e/data/Pseudomonas

# Group by order
cat ASSEMBLY/Pseudomonas.assembly.pass.csv |
    sed -e '1d' |
    tsv-select -d, -f 3 |
    tsv-uniq |
    nwr append stdin -r order |
    tsv-select -f 2 |
    tsv-uniq \
    > order.lst

cat order.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        n_species=$(cat ASSEMBLY/Pseudomonas.assembly.pass.csv |
            sed "1d" |
            tsv-select -d, -f 3 |
            nwr append stdin -r order -r species |
            grep {} |
            tsv-select -f 1,3 |
            tsv-uniq |
            wc -l)

        n_strains=$(cat ASSEMBLY/Pseudomonas.assembly.pass.csv |
            sed "1d" |
            tsv-select -d, -f 3 |
            nwr append stdin -r order |
            grep {} |
            wc -l)

        printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    ' |
    nwr append stdin --id |
    tsv-select -f 5,4,2,3 |
    tsv-sort -k2,2 |
    (echo -e '#tax_id\torder\t#species\t#strains' && cat) |
    mlr --itsv --omd cat
# nwr append will supply tax_id of every order
```

| #tax_id | order                 | #species | #strains |
|---------|-----------------------|----------|----------|
| 1692040 | Acidiferrobacterales  | 3        | 3        |
| 135624  | Aeromonadales         | 18       | 18       |
| 135622  | Alteromonadales       | 52       | 118      |
| 1385    | Bacillales            | 3        | 3        |
| 213849  | Campylobacterales     | 1        | 1        |
| 135615  | Cardiobacteriales     | 2        | 2        |
| 204458  | Caulobacterales       | 1        | 1        |
| 1706369 | Cellvibrionales       | 2        | 4        |
| 51291   | Chlamydiales          | 1        | 1        |
| 135613  | Chromatiales          | 32       | 32       |
| 85007   | Corynebacteriales     | 1        | 1        |
| 91347   | Enterobacterales      | 179      | 179      |
| 1934945 | Immundisolibacterales | 1        | 1        |
| 118969  | Legionellales         | 18       | 18       |
| 135618  | Methylococcales       | 14       | 14       |
| 2887326 | Moraxellales          | 50       | 523      |
| 1775403 | Nevskiales            | 1        | 1        |
| 135619  | Oceanospirillales     | 9        | 18       |
| 1240482 | Orbales               | 3        | 3        |
| 135625  | Pasteurellales        | 34       | 34       |
| 72274   | Pseudomonadales       | 183      | 846      |
| 72273   | Thiotrichales         | 32       | 32       |
| 135623  | Vibrionales           | 46       | 46       |
| 135614  | Xanthomonadales       | 53       | 53       |

- Genus

```bash
cd /mnt/e/data/Pseudomonas

# Group by genus
cat ASSEMBLY/Pseudomonas.assembly.pass.csv |
    sed -e '1d' |
    tsv-select -d, -f 3 |
    tsv-uniq |
    nwr append stdin -r genus |
    tsv-select -f 2 |
    tsv-uniq \
    > genus.lst

cat genus.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        n_species=$(cat ASSEMBLY/Pseudomonas.assembly.pass.csv |
            sed "1d" |
            tsv-select -d, -f 3 |
            nwr append stdin -r genus -r species |
            grep {} |
            tsv-select -f 1,3 |
            tsv-uniq |
            wc -l)

        n_strains=$(cat ASSEMBLY/Pseudomonas.assembly.pass.csv |
            sed "1d" |
            tsv-select -d, -f 3 |
            nwr append stdin -r genus |
            grep {} |
            wc -l)

        printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    ' |
    nwr append stdin --id |
    tsv-select -f 5,4,2,3 |
    tsv-sort -k2,2 |
    tsv-filter --ge 4:10 |
    (echo -e '#tax_id\tgenus\t#species\t#strains' && cat) |
    mlr --itsv --omd cat
# only show the genus contained 10+ strains
# tsv-sort -k2,2 | wc -l
#  171
# totally 171 genera here
```

| #tax_id | genus             | #species | #strains |
|---------|-------------------|----------|----------|
| 469     | Acinetobacter     | 43       | 486      |
| 642     | Aeromonas         | 12       | 12       |
| 226     | Alteromonas       | 16       | 31       |
| 544     | Citrobacter       | 14       | 14       |
| 547     | Enterobacter      | 10       | 10       |
| 262     | Francisella       | 11       | 11       |
| 2745    | Halomonas         | 5        | 13       |
| 445     | Legionella        | 14       | 14       |
| 68      | Lysobacter        | 12       | 12       |
| 475     | Moraxella         | 5        | 35       |
| 122277  | Pectobacterium    | 12       | 12       |
| 53246   | Pseudoalteromonas | 18       | 33       |
| 286     | Pseudomonas       | 172      | 827      |
| 613     | Serratia          | 14       | 14       |
| 22      | Shewanella        | 16       | 52       |
| 662     | Vibrio            | 38       | 38       |
| 338     | Xanthomonas       | 18       | 18       |
| 629     | Yersinia          | 12       | 12       |

- strains

```bash
cd /mnt/e/data/Pseudomonas

# list strains
mkdir -p taxon

rm taxon/* strains.lst *.tmp

# this step will save all strains name into orders
cat ASSEMBLY/Pseudomonas.assembly.pass.csv |
    sed -e '1d' |
    tr "," "\t" |
    tsv-select -f 1,2,3 |
    nwr append stdin -c 3 -r species -r genus -r family -r order |
    parallel --col-sep "\t" --no-run-if-empty --linebuffer -k -j 4 '
        echo {1} >> strains.lst

        echo {4} >> species.tmp
        echo {5} >> genus.tmp
        echo {6} >> family.tmp

        echo {7} >> order.tmp
        echo {1} >> taxon/{7}

        printf "%s\t%s\t%d\t%s\t%s\t%s\t%s\n" {1} {2} {3} {4} {5} {6} {7}
    ' \
    > strains.taxon.tsv
# nwr append [OPTIONS] <infiles> ...
# <infiles> ...: Input filename. [stdin] for standard input
# -c/--column <column>: The column where the IDs are located [default: 1]

cat species.tmp | tsv-uniq > species.lst
cat genus.tmp | tsv-uniq > genus.lst
cat family.tmp | tsv-uniq > family.lst
cat order.tmp | tsv-uniq > order.lst

# Omit strains without protein annotations
for STRAIN in $(cat strains.lst); do
    if ! compgen -G "ASSEMBLY/${STRAIN}/*_protein.faa.gz" > /dev/null; then
        echo ${STRAIN}
    fi
    if ! compgen -G "ASSEMBLY/${STRAIN}/*_cds_from_genomic.fna.gz" > /dev/null; then
        echo ${STRAIN}
    fi
done |
    tsv-uniq \
    > omit.lst

# An example for compgen -G using
compgen -G "ASSEMBLY/Acidif_thio_GCF_001705075_2/*_protein.faa.gz"
#ASSEMBLY/Acidif_thio_GCF_001705075_2/GCF_001705075.2_ASM170507v2_protein.faa.gz

cat omit.lst | wc -l
#  0
# All OK

rm *.tmp
```

> **Explanation**:
> 
> - `compgen`
> 
> ```txt
> Display possible completions depending on the options.
> 
> Intended to be used from within a shell function generating possible completions. If the optional WORD argument is supplied, matches against WORD are generated.
> 
> Exit Status:
> Returns success unless an invalid option is supplied or an error occurs.
> 
> [-G globpat]: Glob pattern 
> ```
> 
> Glob patterns specify sets of filenames with wildcard characters. Unix Bash shell command `mv *.txt textfiles` moves all files with names ending in `.txt` from the current dir to the `textfiles` dir.
> 
> Because we want to find if there is any file of protein seq or CDS seq from ASSEMBLY/${STRAIN}/, so we use `*_protein.faa.gz`.
> 
> `*` is a wildcard character, so use `compgen -G` could identify all possible files with `_protein.faa.gz`.
> 
> For the example showed above, the output is the file dirname/basename if exists.
> 
> The output of `compgen -G` is unnecessary, so `> /dev/null`
> 
> And `!` used to take a test result or exit status opposite

## NCBI taxonomy

> Done by `bp_taxonomy2tree.pl` from BioPerl
>
> - `bp_taxonomy2tree.pl`
> 
> ```txt
> NAME
>   bp_taxonomy2tree - Building a taxonomic tree based on the full lineages
>   of a set of species names
>
> DESCRIPTION
>   This scripts looks up the provided species names in the NCBI Taxonomy
>   database, retrieves their full lineage and puts them in a Newick
>   taxonomic tree displayed on screen.
>
>     bp_taxonomy2tree.pl -s Orangutan -s Gorilla -s Chimpanzee -s Human
>     bp_taxonomy2tree.pl -s Orangutan -s Gorilla -s Chimpanzee -s "Homo Sapiens"
>
> OPTIONS
> -e: use the web-based Entrez tasxonomy database if you do not have the NCBI flatfiles installed
>
> -o: specify your
>
> -a: for the NCBI names file
> ```
>
> - `nw_display`
>
> ```txt
> Displays a tree as a graph, as text or SVG.
>
> Synopsis:
> nw_display [options] <tree filename|->
>
> OPTIONS:
> -: the tree is read on stdin
>
> -b <string>: CSS for branch length lables. [only SVG]
>
> -s: output graph as SVG (default: ASCII graphics). All output is on stdout
>
> -w <number>: graph should be nowider than <number>, measured in characters
>              for text and pixels for SVG. Defaults: 80 (text), 300 (SVG)
>
> -v <number>: number of pixels between leaves
> ```
>
> - `rsvg-convert`
>
> ```txt
> Covert SVG files to other image formats
>
> USAGE:
>   rsvg-convert [FLAGS] [OPTIONS] [FILE]...
> ```

The goal of this part was is to get the taxonomy tree from NCBI.

This tree will supply you the taxonomy info of each species. The following tree can compare to this taxonomy tree and check those anomalously labeled species.

```bash
mkdir -p /mnt/e/data/Pseudomonas/tree
cd /mnt/e/data/Pseudomonas/tree

bp_taxonomy2tree.pl -e \
    $(
        cat ../genus.lst |
            tr " " "_" |
            parallel echo '-s {}'
    ) \
    > ncbi.nwk

nw_display -s -b 'visibility:hidden' -w 600 -v 30 ncbi.nwk |
    rsvg-convert -o Pseudomonas.ncbi.png
```

## Raw phylogenetic tree by MinHash

Using genome to build a tree. It is pretty common that strains from the same species vary widely at the genomic level. The main reason is that the taxa of different strains are seperated by their morphological characteristics, which can hardly reflect their genomic features. So this step could help us find those abnormal strains from genomic level.

### MinHash grouping all strains

MinHash was used and introduced in [plasmid.md](plasmid.md).

- Using MinHash for grouping strains

```bash
mkdir -p /mnt/e/data/Pseudomonas/mash
cd /mnt/e/data/Pseudomonas/mash

for strain in $(cat ../strains.lst ); do
    2>&1 echo "==> ${strain}"

    if [[ -e ${strain}.msh ]]; then
        continue
    fi

    find ../ASSEMBLY/${strain} -name "*_genomic.fna.gz" |
        grep -v "_from_" |
        xargs cat |
        mash sketch -k 21 -s 100000 -p 8 - -I "${strain}" -o ${strain}
done
# mash sketch 
# -k <int>: k-mer 21 
# sketch size: each 100000 will have at most this many non-redundant min-hashes
# -I <path>: ID field for sketch of reads (instead of first sequence ID)

mash triangle -E -p 8 -l <(
    cat ../strains.lst | parallel echo "{}.msh"
    ) \
    > dist.tsv
# mash triangle [options] <seq1> [<seq2>] ...
# Estimate the distance of each input seq to every other input seq.
# Outputs a lower-triangular distance matrix in relaxed Phylip format.
# input files can be fasta or fastq, gzipped or not, or .msh
# -E: output edge list instead of Phylip matrix
# -l: list input, lines in each <query> specify paths to seq files, one per line

ls *.msh | wc -l
#  1952

cat dist.tsv | wc -l
#  1904176
#1952*1951/2*1, C(1952,2)

# fill matrix with lower triangle
tsv-select -f 1-3 dist.tsv |
    (tsv-select -f 2,1,3 dist.tsv && cat) |
    (
        cut -f 1 dist.tsv |
            tsv-uniq |
            parallel -j 1 --keep-order 'echo -e "{}\t{}\t0"' &&
        cat
    ) \
    > dist_full.tsv

tsv-select -f 1-3 dist.tsv | (tsv-select -f 2,1,3 dist.tsv && cat) | wc -l
#  3808352
# actually swap col1 and col2 and then using cat to combine them

cut -f 1 dist.tsv | tsv-uniq | parallel -j 1 --keep-order 'echo -e "{}\t{}\t0"' | wc -l
# other 1951 strains dist to themselves - 0

cat dist_full.tsv | wc -l
#  3810303
# 3808352+1951

cat dist_full.tsv |
    Rscript -e '
        library(readr);
        library(tidyr);
        library(ape);
        pair_dist <- read_tsv(file("stdin"), col_names=F);
        tmp <- pair_dist %>%
            pivot_wider( names_from = X2, values_from = X3, values_fill = list(X3 = 1.0) )
        tmp <- as.matrix(tmp)
        mat <- tmp[,-1]
        rownames(mat) <- tmp[,1]

        dist_mat <- as.dist(mat)
        clusters <- hclust(dist_mat, method = "ward.D2")
        tree <- as.phylo(clusters)
        write.tree(phy=tree, file="tree.nwk")

        group <- cutree(clusters, h=0.4) # k=5
        groups <- as.data.frame(group)
        groups$ids <- rownames(groups)
        rownames(groups) <- NULL
        groups <- groups[order(groups$group), ]
        write_tsv(groups, "groups.tsv")
    '
```

### Tweak the mash tree

> - `nw_reroot`
> 
> ```txt
> (Re)roots a tree on a specified outgroup
> 
> Synopsis:
> nw_reroot [-dhl] <newick trees filename|-> [lable*]
> 
> Input:
> First argument is the name of a file that contains Newick trees, or '-' from stdin
> 
> Further arguments are node lables. If there is at least one label, the tree will be
> re-rooted on their LCA (last common ancestor). If there is no label, the tree is re-
> rooted on the longest branch.
> ```
> 
> - `nw_order`
> 
> ```txt
> Orders nodes according to various criteria, preserving topology
> 
> Synopsis:
> nw_order [-c:hn] <newick trees filename|->
> 
> Input:
> Argument is the name of a file that contains Newick trees, or '-' from stdin
> 
> Output:
> Orders the tree and prints it out on standard output. By default, the ordering
> field is the node's label for leaves, or the first child's order field for inner
> nodes. The tree's topology is not altered: the biological info contained in the
> tree is left intact.
> 
> Options:
> -c <criterion>: specify order criterion. Possible criteria are:
> 'a' (alphanumeric order of labels)
> 'n' (number of descendants: nodes with fewer descendans appear first)
> 'd' (de-ladderize: alternately put nodes with fewer descendants before or after
> those with more)
> ```

```bash
mkdir -p /mnt/e/data/Pseudomonas/tree
cd /mnt/e/data/Pseudomonas/tree

nw_reroot ../mash/tree.nwk Bac_subti_subtilis_168 Sta_aure_aureus_NCTC_8325 |
    nw_order -c n - \
    > mash.reroot.newick
# Two strains are from Firmicutes, and they are both G-positive

# rank::col
ARRAY=(
    'order::7'
    'family::6'
    'genus::5'
    'species::4'
)

rm mash.condensed.map
CUR_TREE=mash.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/withncbi/taxon/condense_tree.sh ${CUR_TREE} ../strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick mash.${GROUP_NAME}.newick
    cat condense.map >> mash.condensed.map

    CUR_TREE=mash.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 600 -v 30 mash.species.newick |
    rsvg-convert -o Pseudomonas.mash.png
```

## Pangenome

> - ppanggolin
>
> Using Partitioned PanGenome Graph Of Linked Neighbors (PPanGGOLiN).
>
> PPanGGOLiN is a software suite used to create and manipulate prokaryotic pangenomes from a set of either genomic DNA sequences or provided genome annotations.
> 
> PPanGGOLiN builds pangenomes through a graphical model and a statistical method to partition gene families in persistent, shell and cloud genomes.
>
> ```txt
> ppanggolin workflow --fasta/--anno ORGANISMS_FASTA_LIST/ORGANISMS_ANNOTATION_LIST
> It uses parameters that we found to be generally the best when working with species pangenomes.
>
> The file ORGANISMS_FASTA_LIST is a tsv-separated file with the following organisation:
>   1. The first col contains a unique organism name (without whitespace)
>   2. The second col contains the path to the associated FASTA/GFF file
>   3. (Only in fasta) Circular contig identifiers are indicated
>   4. Each line represents an organism
> ```

```bash
cd /mnt/e/data/Pseudomonas

mkdir -p pangenome

cat strains.lst |
    grep "^Pseudom_aeru_" |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        GBFF=$(compgen -G "ASSEMBLY/{}/*_genomic.gbff.gz")

        if [ "${GBFF}" == "" ]; then
            exit;
        fi

        echo {}
        compgen -G "ASSEMBLY/{}/*_genomic.gbff.gz"
    ' |
    paste - - \
    > pangenome/Pseudom_aeru.gbff.list
# grep "^Pseudom_aeru_": will get all P. aeruginosa
# after paste, it will be the name, protein file path

wc -l < pangenome/Pseudom_aeru.gbff.list
#  391

ppanggolin workflow --anno pangenome/Pseudom_aeru_1.gbff.list --cpu 8 -o pangenome/Pseudom_aeru
# ppanggolin workflow
# .gbff.list example could be seen above
```

## Collect proteins

> **The goal of this part**:
> 
> - First, the following steps are focused on proteins among all strains. So do some rough preliminary statistics on them.
> - Second, merge all proteins to a file could reduce the complexity of the following steps in manipulating different kinds of proteins.


- `all.pro.fa`

Merge all proteins to a file and using accession to get uniq one.

```bash
cd /mnt/e/data/Pseudomonas

mkdir -p PROTEINS

find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
    sort |
    grep 'ASSEMBLY/' |
    wc -l
#  1958

find ASSEMBLY -type f -name "*_protein.faa.gz" |
    wc -l
#  1958

cat strains.lst |
    wc -l
#  1952
# All strains after filtering

for STRAIN in $(cat strains.lst); do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz
done \
    > PROTEINS/all.pro.fa

cat PROTEINS/all.pro.fa |
    perl -nl -e '
        BEGIN { our %seen; our $h; }

        if (/^>/) {
            $h = (split(" ", $_))[0];
            $seen{$h}++;
            $_ = $h;
        }
        print if $seen{$h} == 1;
    ' \
    > PROTEINS/all.uniq.fa
# $h in perl represented protein accessions right after '>'
# $seen{$h} == 1 will gives out those proteins when they first occurred
# occurrence will be recorded by hash %seen

# counting proteins
cat PROTEINS/all.pro.fa |
    grep "^>" |
    wc -l
#  8754303

cat PROTEINS/all.pro.fa |
    grep "^>" |
    tsv-uniq |
    wc -l
#  3995205

# annotations may be different
cat PROTEINS/all.uniq.fa |
    grep "^>" |
    wc -l
#  3944568
# which means the same accession annotated differently 

# ribonuclease was an example about how to grep proteins needed from the whole protein fasta
cat PROTEINS/all.pro.fa |
    grep "ribonuclease" |
    grep -v "deoxyribonuclease" |
    perl -nl -e 's/^>\w+\.\d+\s+//g; print' | # removed accession
    perl -nl -e 's/\s+\[.+?\]$//g; print' | # removed [SPECIES]
    perl -nl -e 's/MULTISPECIES: //g; print' |
    sort |
    uniq -c |
    sort -nr
 # only grep RNase, so the DNase should be discarded
```

- `all.replace.fa`

Replace all proteins' name into strains_accession format. 

```bash
cd /mnt/e/data/Pseudomonas

rm PROTEINS/all.strain.tsv
for STRAIN in $(cat strains.lst); do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
        grep "^>" |
        cut -d" " -f 1 |
        sed "s/^>//" |
        STRAIN=${STRAIN} perl -nl -e '
            $n = $_;
            $s = $n;
            $s =~ s/\.\d+//;
            printf qq{%s\t%s_%s\t%s\n}, $n, $ENV{STRAIN}, $s, $ENV{STRAIN};
        ' \
    > PROTEINS/${STRAIN}.replace.tsv

    cut -f 2,3 PROTEINS/${STRAIN}.replace.tsv >> PROTEINS/all.strain.tsv

    faops replace -s ASSEMBLY/${STRAIN}/*_protein.faa.gz <(cut -f 1,2 PROTEINS/${STRAIN}.replace.tsv) stdout

    rm PROTEINS/${STRAIN}.replace.tsv
done \
    > PROTEINS/all.replace.fa
# In perl, environment variables in a special hash named %ENV
# So you could use $ENV{STRAIN} to access STRAIN in ENV list
# faops replace: Replace headers from a FA file
# -s: ouly output seqs in the list, like `faops some`
# <replace.tsv> is a tab-separated file containing two fields
# original_name \t replace_name
# After running commands above, all proteins' name were changed from accessions to names

cat PROTEINS/all.replace.fa |
    grep "^>" |
    wc -l
#  8754303
# the same to all.pro.fa, so all protein names were replaced by strain_accession(pro)

(echo -e "#name\tstrain" && cat PROTEINS/all.strain.tsv)  \
    > temp &&
    mv temp PROTEINS/all.strain.tsv

faops size PROTEINS/all.replace.fa > PROTEINS/all.replace.sizes

(echo -e "#name\tsize" && cat PROTEINS/all.replace.sizes) > PROTEINS/all.size.tsv

rm PROTEINS/all.replace.sizes
```

- `all.info.tsv`

Combine all info into one file - including name, strain, size, and annotation.

```bash
cd /mnt/e/data/Pseudomonas

for STRAIN in $(cat strains.lst); do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
        grep "^>" |
        sed "s/^>//" |
        perl -nl -e '/\[.+\[/ and s/\[/\(/; print' |
        perl -nl -e '/\].+\]/ and s/\]/\)/; print' |
        perl -nl -e 's/\s+\[.+?\]$//g; print' |
        perl -nl -e 's/MULTISPECIES: //g; print' |
        STRAIN=${STRAIN} perl -nl -e '
            /^(\w+)\.\d+\s+(.+)$/ or next;
            printf qq{%s_%s\t%s\n}, $ENV{STRAIN}, $1, $2;
        '
done \
    > PROTEINS/all.annotation.tsv

cat PROTEINS/all.annotation.tsv |
    wc -l
#  8754303

(echo -e "#name\tannotation" && cat PROTEINS/all.annotation.tsv) \
    > temp &&
    mv temp PROTEINS/all.annotation.tsv

# check differences
cat PROTEINS/all.size.tsv |
    grep -F -f <(cut -f 1 PROTEINS/all.annotation.tsv) -v
# -F/--fixed-strings: comparing strings directly, whole strings as PATTERNS
# -f/--file=FILE: take PATTERNS from FILE
# nothing showed on screen, so col1 from both files were the same 

tsv-join \
    PROTEINS/all.strain.tsv \
    --data-fields 1 \
    -f PROTEINS/all.size.tsv \
    --key-fields 1 \
    --append-fields 2 \
    > PROTEINS/all.strain_size.tsv

tsv-join \
    PROTEINS/all.strain_size.tsv \
    --data-fields 1 \
    -f PROTEINS/all.annotation.tsv \
    --key-fields 1 \
    --append-fields 2 \
    > PROTEINS/all.info.tsv

cat PROTEINS/all.info.tsv |
    wc -l
#  8754304
```

## Phylogenetics with 40 single-copy genes

This part is using proteins inside 40 single-copy genes to generate a protein tree.

### Find correspoding proteins by `hmmsearch`

- Download HMM models as described in `hmm/README.md`
- The `E_VALUE` was manually adjusted to 1e-20 to reach a balance between sensitivity and speciality.

HMM ([Hidden Markov Model](https://en.wikipedia.org/wiki/Hidden_Markov_model)) is a statistical Markov model in which the system being modeled is assumed to be a Markov process - call it $X$ - with unobservable ("hidden") states.

As part of the definition, HMM requires that there be an observable process $Y$ whose outcomes are "influenced" by the outcomes of $X$ in a known way.

Different transduction matrix would give out different multi-alignment results.

> - hmmsearch
> 
> ```txt
> hmmsearch :: search profile(s) against a sequence database
> HMMER 3.3.2 (Nov 2020); https://hmmer.org/
> 
> Usage: hmmsearch [options] <hmmfile> <seqdb>
> 
> Either the query <hmmfile> or the target <seqdb> may be '-', in which case the query
> profile or target database input will be read from a <stdin> pipe instead of from a file.
> 
> Options controlling reporting thresholds:
> -E <x>:  report sequences <= this E-value threshold in output  [10.0]  (x>0)
> --domE <x>: report domains <= this E-value threshold in output  [10.0]  (x>0)
> 
> Options directing output:
> --noali: don't output alignments, so output is smaller
> --notextw: unlimit ASCII text output line width
> ```
> 
> `hmmsearch` may take minutes or even hours to run, dpending on thee size of the seq database.
> 
> The output consists of four sections: a ranked list of the best scoring sequences, a ranked list of the best scoring domains, alignments for all the best scoring domains, and a histogram of the scores.

- An example

```bash
gzip -dcf ASSEMBLY/Acidif_thio_GCF_001705075_2/*_protein.faa.gz |
    hmmsearch --noali --notextw ~/data/HMM/scg40/bacteria_and_archaea_dir/BA00004.hmm - |
    grep '>>' |
    perl -nl -e '/>>\s+(\S+)/ and print $1'
#WP_083995962.1
#WP_065970328.1
#WP_065969050.1
#WP_065969049.1
#WP_065972062.1
#WP_065971275.1
#WP_228579725.1
#WP_065971817.1

gzip -dcf ASSEMBLY/Acidif_thio_GCF_001705075_2/*_protein.faa.gz |
    hmmsearch -E 1e-20 --domE 1e-20 --noali --notextw ~/data/HMM/scg40/bacteria_and_archaea_dir/BA00004.hmm - |
    grep '>>' |
    perl -nl -e '/>>\s+(\S+)/ and print $1'
#WP_083995962.1
#WP_065970328.1
#WP_065969050.1
#WP_065969049.1
#WP_065972062.1
#WP_065971275.1

# the example showed above was extracted domain accession from the hmmsearch results
# meanwhile, strict E-value could filter and then leave less proteins more close to markers.  
```

```bash
E_VALUE=1e-20

cd /mnt/e/data/Pseudomonas

# Find all genes
for marker in BA000{01..40}; do
    >&2 echo "==> marker [${marker}]"

    mkdir -p PROTEINS/${marker}

    for ORDER in $(cat order.lst); do
        >&2 echo "==> ORDER [${ORDER}]"

        cat taxon/${ORDER} |
            parallel --no-run-if-empty --linebuffer -k -j 8 "
                gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                    hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw ~/data/HMM/scg40/bacteria_and_archaea_dir/${marker}.hmm - |
                    grep '>>' |
                    perl -nl -e ' m{>>\s+(\S+)} and printf qq{%s\t%s\n}, \$1, {}; '
            " \
            > PROTEINS/${marker}/${ORDER}.replace.tsv
    done

    echo
done
# will get accessions of protein and strains from an order into a replace.tsv
```

### Create a valid marker gene list

```bash
for marker in BA000{01..40}; do
    >&2 echo "==> marker [${marker}]"
    
    for ORDER in $(cat order.lst); do
        if [ -s PROTEINS/${marker}/${ORDER}.replace.tsv ]; then
            copy=$(cat PROTEINS/${marker}/${ORDER}.replace.tsv |
                      tsv-summarize --g 2 --count |
                      tsv-summarize --g 2 --count |
                      sort -r -nk 2,2 |
                      head -n 1 |
                      tsv-select -f 1
                  )

            echo "==> ORDER [${ORDER}] copy=${copy}"

            cat PROTEINS/${marker}/${ORDER}.replace.tsv |
                tsv-summarize -g 2 --count |
                tsv-filter --ne 2:${copy}
        else
            echo "==> ORDER [${ORDER}] no copy"
        fi
    done

    echo
done
```

After running the command above, those marker genes are indeed unstable among species. We remove them afterwards.

- `hmmsearch` may identify more than one copy for some marker genes
  - BA00004: translation initiation factor EF-2
  - BA00005: translation initiation factor IF-2
  - BA00008: signal recognition particle protein
- Acinetobacter
  - BA00028

Compare proteins and strains

```bash
cd /mnt/e/data/Pseudomonas

for marker in BA000{01..03} BA000{06..07} BA000{09..27} BA000{29..40}; do
    echo ${marker}
done > marker.lst

for marker in $(cat marker.lst); do
    echo "==> marker [${marker}]"

    for ORDER in $(cat order.lst); do
        cat PROTEINS/${marker}/${ORDER}.replace.tsv |
            cut -f 2 |
            diff - taxon/${ORDER}
    done

    echo
done

# Part of the output:
#==> marker [BA00038]
#434d433
#< Acin_pittii_PHEA_2
#0a1
#> Act_del_GCF_900638385_1
#13a15
#> Hae_aeg_GCF_900475885_1
#31a34
#> Ve_pul_GCF_013377275_1
#2a3
#> Are_dae_GCF_007993735_1
#1a2
#> Alt_aus_GCF_000934525_1
#4a5
#> Leg_antarctica_GCF_011764505_1
#6a8,9
#> Leg_gee_GCF_004571195_1
#> Leg_hac_GCF_000953655_1
#7a11
#> Leg_jor_GCF_900637635_1
#699a700
#> Pseudom_puti_GCF_014854555_1
#0a1
#> Chl_tracho_D_UW_3_CX

# diff: Compare FILES line by line
# diff file1 file2
# diff output (so called 'normal diff'):
# < denotes lines in file1 (in this case means -, that is order replace by hmmsearch)
# > denotes lines in file2 (in this case means taxon/ files with the corresponding name)
# d stands for deletion, a stands for adding, c stands for changing
# 434d433 means the 434 line number in file1 was deleted and has the line number 433 in file2
# 4a5 means starting form line number 4 in file1 added and this added line is the number 5 in file2
```

### Align and concat marker genes to create species tree

- Strains within a species share a large proportion of identical protein sequences.
- Use `trimal -automated1` to remove gaps.

> - `trimAL`
> 
> ```txt
> trimAl is a tool for the automated removal of spurious sequences or poorly aligned regions
> from a multiple sequence alignment.
> 
> Basic usage:
> trimal -in <inputfile> -out <outputfile> -(other options)
> 
> Options:
> -automated1: Use a heuristic selection of the automatic method based on similarity statistics.
> (Optimized for Maximum Likelihood phylogenetic tree reconstruction).
> ```
> 
> - `muscle`
> 
> ```txt
> MUSCLE stands for MUltiple Sequence Comparison by Log-Expectation.
> 
> Basic usage:
> muscle -in <inputfile> -out <outputfile>
> 
> Options:
> -quiet: Do not write progress messages to stderr
> ```
> 
> - `FastTree`
> 
> ```txt
> FastTree infers approximately-maximum-likelihood phylogenetic trees
> from alignments of nucleotide or protein sequences.
> 
> Basic usage:
> FastTree protein_alignment > tree
> FastTree -nt nucleotide_alignment > tree
> 
> Options:
> --fastest: speed up the neighbor foining phase & reduce memory usage
> --noml: to turn off maximum-likelihood
> ```

```bash
cd /mnt/e/data/Pseudomonas

# Extract sequences
cat marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo "==> marker [{}]"

        for ORDER in $(cat order.lst); do
            cat PROTEINS/{}/${ORDER}.replace.tsv
        done \
            > PROTEINS/{}/{}.replace.tsv

        faops some PROTEINS/all.uniq.fa <(
            cat PROTEINS/{}/{}.replace.tsv |
                cut -f 1 |
                tsv-uniq
            ) stdout \
            > PROTEINS/{}/{}.pro.fa
    '
# cat a marker sub-dir into one file, and then using faops some to extract protein fasta

# Align each markers with `muscle`
cat marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo "==> marker [{}]"

        muscle -quiet -in PROTEINS/{}/{}.pro.fa -out PROTEINS/{}/{}.aln.fa
    '

for marker in $(cat marker.lst); do
    >&2 echo "==> marker [${marker}]"

    # 1 name to many names
    cat PROTEINS/${marker}/${marker}.replace.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 4 "
            faops replace -s PROTEINS/${marker}/${marker}.aln.fa <(echo {}) stdout
        " \
        > PROTEINS/${marker}/${marker}.replace.fa
done
# because .replace.tsv got from the above command lines contained two cols
# so faops replace will replace the accession to strain_pro

# Concat marker genes
for marker in $(cat marker.lst); do
    # sequences in one line
    faops filter -l 0 PROTEINS/${marker}/${marker}.replace.fa stdout

    # empty line for .fas
    echo
done \
    > PROTEINS/scg40.aln.fas

fasops concat PROTEINS/scg40.aln.fas strains.lst -o PROTEINS/scg40.aln.fa
# fasops concat: Concatenate sequence pieces in blocked fasta files.
# so actually sequences from 40 maker genes will be concatenated into one 

# Trim poorly aligned regions with `TrimAl`
trimal -in PROTEINS/scg40.aln.fa -out PROTEINS/scg40.trim.fa -automated1

# FastTree produces NJ trees to simulate ML ones
FastTree PROTEINS/scg40.trim.fa > PROTEINS/scg40.trim.newick
```

### Tweak the concat tree

```bash
cd /mnt/e/data/Pseudomonas/tree

nw_reroot ../PROTEINS/scg40.trim.newick Bac_subti_subtilis_168 Sta_aure_aureus_NCTC_8325 |
    nw_order -c n - \
    > scg40.reroot.newick

# rank::col
ARRAY=(
#    'order::7'
#    'family::6'
    'genus::5'
    'species::4'
)

rm scg40.condensed.map
CUR_TREE=scg40.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/withncbi/taxon/condense_tree.sh ${CUR_TREE} ../strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick scg40.${GROUP_NAME}.newick
    cat condense.map >> scg40.condensed.map

    CUR_TREE=scg40.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 600 -v 30 scg40.species.newick |
    rsvg-convert -o Pseudomonas.scg40.png
```

## Phlogenetics with bac120

Bac120 was introduced in [hmm.md](hmm.md). This bac120 contained whole length proteins from TIGRFAM database and some domains. Now, this method has gradually occupied the mainstream in the field of generating phylogenetic tree from proteins. The following steps are relied more on results from bac120.

### Find corresponding proteins by `hmmsearch`

```bash
E_VALUE=1e-20

cd /mnt/e/data/Pseudomonas

# Find all genes
for marker in $(cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1); do
    >&2 echo "==> marker [${marker}]"

    mkdir -p PROTEINS/${marker}

    for ORDER in $(cat order.lst); do
        >&2 echo "==> ORDER [${ORDER}]"

        cat taxon/${ORDER} |
            parallel --no-run-if-empty --linebuffer -k -j 8 "
                gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                    hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw ~/data/HMM/bac120/HMM/${marker}.HMM - |
                    grep '>>' |
                    perl -nl -e ' m{>>\s+(\S+)} and printf qq{%s\t%s\n}, \$1, {}; '
            " \
            > PROTEINS/${marker}/${ORDER}.replace.tsv
    done

    echo
done
```

### Align and concat marker genes to create species tree

```bash
cd /mnt/e/data/Pseudomonas

cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        for ORDER in $(cat order.lst); do
            cat PROTEINS/{}/${ORDER}.replace.tsv
        done |
            wc -l
    ' |
    tsv-summarize --quantile 1:0.25,0.5,0.75
#  1948    1951    3066

cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo {}
        for ORDER in $(cat order.lst); do
            cat PROTEINS/{}/${ORDER}.replace.tsv
        done |
            wc -l
    ' |
    paste - - |
    tsv-filter --invert --ge 2:1500 --le 2:2500 |
    cut -f 1 \
    > PROTEINS/bac120.omit.lst

# Extract sequences
# Multiple copies slow down the alignment process
cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    grep -v -Fx -f PROTEINS/bac120.omit.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo "==> marker [{}]"

        for ORDER in $(cat order.lst); do
            cat PROTEINS/{}/${ORDER}.replace.tsv
        done \
            > PROTEINS/{}/{}.replace.tsv

        faops some PROTEINS/all.uniq.fa <(
            cat PROTEINS/{}/{}.replace.tsv |
                cut -f 1 |
                tsv-uniq
            ) stdout \
            > PROTEINS/{}/{}.pro.fa
    '

# Align each markers with muscle
cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        >&2 echo "==> marker [{}]"
        if [ ! -s PROTEINS/{}/{}.pro.fa ]; then
            exit
        fi
        if [ -s PROTEINS/{}/{}.aln.fa ]; then
            exit
        fi

        muscle -quiet -in PROTEINS/{}/{}.pro.fa -out PROTEINS/{}/{}.aln.fa
    '
# if [ -s FILE ]: True if FILE exists and has a size greater than zero.

for marker in $(cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1); do
    >&2 echo "==> marker [${marker}]"
    if [ ! -s PROTEINS/${marker}/${marker}.pro.fa ]; then
        continue
    fi

    # 1 name to many names
    cat PROTEINS/${marker}/${marker}.replace.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 4 "
            faops replace -s PROTEINS/${marker}/${marker}.aln.fa <(echo {}) stdout
        " \
        > PROTEINS/${marker}/${marker}.replace.fa
done

# Concat marker genes
for marker in $(cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1); do
    if [ ! -s PROTEINS/${marker}/${marker}.pro.fa ]; then
        continue
    fi

    # sequences in one line
    faops filter -l 0 PROTEINS/${marker}/${marker}.replace.fa stdout

    # empty line for .fas
    echo
done \
    > PROTEINS/bac120.aln.fas

fasops concat PROTEINS/bac120.aln.fas strains.lst -o PROTEINS/bac120.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in PROTEINS/bac120.aln.fa -out PROTEINS/bac120.trim.fa -automated1

faops size PROTEINS/bac120.*.fa |
    tsv-uniq -f 2 |
    cut -f 2
#52756
#25526

# To make it faster
FastTree -fastest -noml PROTEINS/bac120.trim.fa > PROTEINS/bac120.trim.newick
```

### Tweak the concat tree

```bash
cd /mnt/e/data/Pseudomonas/tree

nw_reroot ../PROTEINS/bac120.trim.newick Bac_subti_subtilis_168 Sta_aure_aureus_NCTC_8325 |
    nw_order -c n - \
    > bac120.reroot.newick

# rank::col
ARRAY=(
#    'order::7'
#    'family::6'
    'genus::5'
    'species::4'
)

rm bac120.condensed.map
CUR_TREE=bac120.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/withncbi/taxon/condense_tree.sh ${CUR_TREE} ../strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick bac120.${GROUP_NAME}.newick
    cat condense.map >> bac120.condensed.map

    CUR_TREE=bac120.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 600 -v 30 bac120.species.newick |
    rsvg-convert -o Pseudomonas.bac120.png
```

After all the different trees, we could compare 3 different trees.

## Protein domains and families

### Proteins in Pseudomonas strains

- [`GO:0005975` carbohydrate metabolic process](https://www.ebi.ac.uk/QuickGO/GTerm?id=GO:0005975)
  - http://www.cazy.org/Glycoside-Hydrolases.html

- [`GO:0016837` carbon-oxygen lyase activity, acting on polysaccharides](https://www.ebi.ac.uk/QuickGO/GTerm?id=GO:0016837)
  - http://www.cazy.org/Polysaccharide-Lyases.html

- Search `https://www.pseudomonas.com/goterms` for `GO:0005975`
  - 107 - Pseudomonas aeruginosa PAO1
  - 7360 - Pseudomonas putida KT2440
  - 479 - Pseudomonas chlororaphis subsp. aureofaciens 30-84
  - 116 - Pseudomonas fluorescens SBW25
  - 113 - Pseudomonas protegens Pf-5
  - 123 - Pseudomonas stutzeri A1501
  - 112 - Pseudomonas syringae pv. syringae B728a
  - 114 - Pseudomonas savastanoi pv. phaseolicola 1448A
  - 117 - Pseudomonas entomophila L48
  - 109 - Pseudomonas aeruginosa UCBPP-PA14

- The `E_VALUE` was adjusted to 1e-5 to capture all possible sequences. Because in this part, the goal is to find any possible domain in all strains, so the E_VALUE should be adapted for this goal. Greater value in E_VALUE will give us more domains.

- Using `GO` to extract their pro_seqs from the gff of ref_strain.

```bash
cd /mnt/e/data/Pseudomonas/
mkdir -p /mnt/e/data/Pseudomonas/DOMAINS/ref_strain

for ID in 107 7360 479 116 113 123 112 114 117 109 ; do
    URL=$(printf 'https://www.pseudomonas.com/goterms/list?accession=GO:0005975&strain_id=%d&format=TAB' $ID)
    curl -L ${URL}
done  |
    sed '1d' |
    tsv-select -f 1 \
    > DOMAINS/ref_strain/locus.lst

# this step will split ID of genes from gff using locus_tags as key strings
gzip -dcf \
    ASSEMBLY/Pseudom_aeru_PAO1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_puti_KT2440_GCF_000007565_2/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_chl_aureofaciens_30_84_GCF_000281915_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_fluo_SBW25_GCF_000009225_2/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_prot_Pf_5_GCF_000012265_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_stu_A1501_GCF_000013785_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_syr_pv_syringae_B728a_GCF_000012245_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_sav_pv_phaseolicola_1448A_GCF_000012205_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_entomophi_L48_GCF_000026105_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_aeru_UCBPP_PA14_GCF_000014625_1/*_genomic.gff.gz |
    grep -v "^#" |
    grep -F -w -f DOMAINS/ref_strain/locus.lst |
    tsv-filter --str-eq 3:gene |
    perl -nl -e 'print $1 if /\bID=(.+?);/i' \
    > DOMAINS/ref_strain/gene.lst
# grep
# -w/--word-regexp: match only whole words
# -F/--fixed-strings: PATTERNS are strings

# this step will split Names of CDS from gff using gene IDs as key strings
gzip -dcf \
    ASSEMBLY/Pseudom_aeru_PAO1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_puti_KT2440_GCF_000007565_2/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_chl_aureofaciens_30_84_GCF_000281915_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_fluo_SBW25_GCF_000009225_2/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_prot_Pf_5_GCF_000012265_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_stu_A1501_GCF_000013785_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_syr_pv_syringae_B728a_GCF_000012245_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_sav_pv_phaseolicola_1448A_GCF_000012205_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_entomophi_L48_GCF_000026105_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_aeru_UCBPP_PA14_GCF_000014625_1/*_genomic.gff.gz |
    grep -v "^#" |
    grep -F -w -f DOMAINS/ref_strain/gene.lst |
    tsv-filter --str-eq 3:CDS |
    perl -nl -e 'print $1 if /\bName=(.+?);/i' \
    > DOMAINS/ref_strain/pro.lst

wc -l DOMAINS/ref_strain/*.lst
#  497 DOMAINS/ref_strain/gene.lst
#  513 DOMAINS/ref_strain/locus.lst
#  497 DOMAINS/ref_strain/pro.lst

faops some PROTEINS/all.uniq.fa DOMAINS/ref_strain/pro.lst stdout \
    > DOMAINS/ref_strain/pro.fa

faops size DOMAINS/ref_strain/pro.fa | wc -l
#  495
# 2 miss from pro.lst

# prepare an HMM database for faster hmmscan searches
# hmmpress ~/data/HMM/PFAM/Pfam-A.hmm

E_VALUE=1e-5

hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw \
    ~/data/HMM/PFAM/Pfam-A.hmm DOMAINS/ref_strain/pro.fa |
    grep '>>' |
    perl -nl -e '/>>\s+(\S+)/ and print $1' |
    tee DOMAINS/ref_strain/domain.lst
# the same arguments with hmmsearch
# --noali: don't output alignments, so output is smaller
# --notextw: unlimit ASCII text output line width
# -E <x>: report models <= this E-value threshold in output  [10.0]  (x>0)
# --domE <x>: report domains <= this E-value threshold in output  [10.0]  (x>0)

# query ID against .hmm.dat to find AC
cat ~/data/HMM/PFAM/Pfam-A.hmm.dat |
    grep -F -w -f DOMAINS/ref_strain/domain.lst -A 2 |
    grep -E ' (ID|AC) ' |
    perl -nl -e 'print substr($_, 10) ' |
    paste -d $'\t' - - |
    perl -nlp -e 's/\.\d+$//g' |
    tsv-select -f 2,1 \
    > DOMAINS/ref_strain/pfam_domain.tsv
# grep -A/--after-context=NUM: print NUM lines of trailing context
# meaning the output contained grep and more 2 lines after

wc -l < DOMAINS/ref_strain/pfam_domain.tsv
#  107
```

The results contain all pfam_domain after scanning by `hmmscan`, and those 107 domains are all related to `GO:0005975` carbohydrate metabolic process. All domains are generated from ref_strain.

### Scrap PFAM domains

- Perform keyword search in `pfam` and save the result page as `html only`.
  - [GO:0005975](https://pfam.xfam.org/search/keyword?query=GO%3A0005975)
  - [Glyco_hyd](http://pfam.xfam.org/search/keyword?query=Glyco_hyd)

> - `pup`
> 
> ```txt
> Usage:
> cat index.html | pup [flags] '[selectors] [optional display function]'
> ```
> 
> `tldr pup`
> 
> ```txt
> Command-line HTML parsing tool.
> More information: <https://github.com/ericchiang/pup>.
> 
> - Print all text from the filtered HTML elements and their children:
>     cat index.html | pup 'div text{}'
> ```

```bash
cd /mnt/e/data/Pseudomonas/

cp DOMAINS/ref_strain/pfam_domain.tsv raw.tsv

cat GO_0005975.htm |
    pup 'table.resultTable tr td text{}' |
    grep '\S' |
    paste -d $'\t' - - - - - |
    tsv-select -f 2-4 \
    >> raw.tsv

cat Glyco_hyd.htm |
    pup 'table.resultTable tr td text{}' |
    grep '\S' |
    paste -d $'\t' - - - - - |
    tsv-select -f 2-4 \
    >> raw.tsv
# \S: is a negated \s; it represents any non-whitespace character, [^\s]

#* [carbohydrate metabolic](https://pfam.xfam.org/search/keyword?query=carbohydrate+metabolic)
#cat carbohydrate_metabolic.htm |
#    pup 'table.resultTable tr td text{}' |
#    grep '\S' |
#    paste -d $'\t' - - - - - - - - |
#    tsv-select -f 2-4 \
#    >> raw.tsv

#* [Glycosyl_hydrolase](https://pfam.xfam.org/search/keyword?query=Glycosyl+hydrolase)
#cat Glycosyl_hydrolase.htm |
#    pup 'table.resultTable tr td text{}' |
#    grep '\S' |
#    paste -d $'\t' - - - - - - - - - |
#    tsv-filter --ne 5:10000000 | # Text fields of Pfam entries are not empty
#    tsv-select -f 2-4 \
#    >> raw.tsv

wc -l < raw.tsv
#  290

cat raw.tsv |
    tsv-filter --str-not-in-fld 2:"DUF" |
    tsv-uniq -f 1 |
    tsv-sort -k2,2 \
    > pfam_domain.tsv
# DUF represents domain with unknown function

wc -l < pfam_domain.tsv
#  216

mkdir -p ~/data/Pseudomonas/DOMAINS/HMM

cat pfam_domain.tsv |
    parallel --col-sep "\t" --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo {2}
        curl -L https://pfam.xfam.org/family/{1}/hmm > DOMAINS/HMM/{2}.hmm
    '

find DOMAINS/HMM -type f -name "*.hmm" |
    wc -l
#  216
```

Basically, using 2 different methods, we got totally 216 domains from both PFAM and Pseudomonas GO term after removed those domains whose function have been unknown yet.

`pfam_domain.tsv` included all domains that related to carbohydrate metabolic process.

### Scan every domain

Main goal is to find those domains from our strains by `hmmsearch`.

- The `E_VALUE` was adjusted to 1e-5 to capture all possible sequences

```bash
E_VALUE=1e-5

cd /mnt/e/data/Pseudomonas/

for domain in $(cat pfam_domain.tsv | cut -f 2 | sort); do
    >&2 echo "==> domain [${domain}]"

    if [ -e DOMAINS/${domain}.replace.tsv ]; then
        continue;
    fi

    for ORDER in $(cat order.lst); do
        >&2 echo "==> ORDER [${ORDER}]"

        cat taxon/${ORDER} |
            parallel --no-run-if-empty --linebuffer -k -j 8 "
                gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                    hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw DOMAINS/HMM/${domain}.hmm - |
                    grep '>>' |
                    perl -nl -e '
                        m{>>\s+(\S+)} or next;
                        \$n = \$1;
                        \$s = \$n;
                        \$s =~ s/\.\d+//;
                        printf qq{%s\t%s_%s\n}, \$n, {}, \$s;
                    '
            "
    done \
        > DOMAINS/${domain}.replace.tsv

    >&2 echo
done

ls DOMAINS/*.replace.tsv | wc -l
#  216
# now each DOMAINS/${domain}.replace.tsv has all proteins info that has the ${domain}
# all 216 domains have their own tsv contained proteins from all strains

for domain in $(cat pfam_domain.tsv | cut -f 2 | sort); do
    wc -l DOMAINS/${domain}.replace.tsv
done |
    datamash reverse -W |
    tsv-filter --ge 2:2000 |
    tsv-sort -k2,2nr |
    (echo -e "Domain\tCount" && cat) |
    mlr --itsv --omd cat
# datamash reverse: reverse cols number, which means col-end -> col-start
# -W/--whitespace: use whitespace (one or more spaces and/or tabs) for field delimiters
# this step could give you each domains count
# a strain could have many proteins contained a very specific domain
```

| Domain                              | Count |
|-------------------------------------|-------|
| DOMAINS/Epimerase.replace.tsv       | 40284 |
| DOMAINS/Hydrolase.replace.tsv       | 35532 |
| DOMAINS/HAD.replace.tsv             | 22725 |
| DOMAINS/CBS.replace.tsv             | 21493 |
| DOMAINS/GDP_Man_Dehyd.replace.tsv   | 21275 |
| DOMAINS/HAD_2.replace.tsv           | 20959 |
| DOMAINS/RmlD_sub_bind.replace.tsv   | 19014 |
| DOMAINS/Glycos_transf_2.replace.tsv | 18461 |
| DOMAINS/F420_oxidored.replace.tsv   | 15900 |
| DOMAINS/Glyco_tranf_2_3.replace.tsv | 15734 |
| DOMAINS/3Beta_HSD.replace.tsv       | 15218 |
| DOMAINS/Hydrolase_3.replace.tsv     | 14970 |
| DOMAINS/Hydrolase_like.replace.tsv  | 14372 |
| DOMAINS/PfkB.replace.tsv            | 11621 |
| DOMAINS/SIS.replace.tsv             | 11589 |
| DOMAINS/AAA_33.replace.tsv          | 10988 |
| DOMAINS/Glyco_trans_4_4.replace.tsv | 10216 |
| DOMAINS/Alpha-amylase.replace.tsv   | 8387  |
| DOMAINS/AP_endonuc_2.replace.tsv    | 7594  |
| DOMAINS/ApbA.replace.tsv            | 7173  |
| DOMAINS/Glyco_trans_2_3.replace.tsv | 7143  |
| DOMAINS/Phos_pyr_kin.replace.tsv    | 6989  |
| DOMAINS/CTP_transf_like.replace.tsv | 6561  |
| DOMAINS/Polysacc_deac_1.replace.tsv | 6150  |
| DOMAINS/FGGY_C.replace.tsv          | 6079  |
| DOMAINS/PGM_PMM_I.replace.tsv       | 5353  |
| DOMAINS/FGGY_N.replace.tsv          | 5288  |
| DOMAINS/PGM_PMM_II.replace.tsv      | 5240  |
| DOMAINS/PGM_PMM_III.replace.tsv     | 5234  |
| DOMAINS/PGM_PMM_IV.replace.tsv      | 4831  |
| DOMAINS/Glyco_transf_21.replace.tsv | 4791  |
| DOMAINS/NAD_Gly3P_dh_N.replace.tsv  | 4291  |
| DOMAINS/SKI.replace.tsv             | 3889  |
| DOMAINS/QRPTase_C.replace.tsv       | 3824  |
| DOMAINS/Glyco_hydro_3.replace.tsv   | 3698  |
| DOMAINS/Aldose_epim.replace.tsv     | 3536  |
| DOMAINS/CBM_48.replace.tsv          | 3407  |
| DOMAINS/DctQ.replace.tsv            | 3160  |
| DOMAINS/LamB_YcsF.replace.tsv       | 2784  |
| DOMAINS/Glyco_transf_28.replace.tsv | 2446  |
| DOMAINS/Glyco_tran_28_C.replace.tsv | 2409  |
| DOMAINS/TAL_FSA.replace.tsv         | 2340  |
| DOMAINS/F_bP_aldolase.replace.tsv   | 2239  |
| DOMAINS/Ribul_P_3_epim.replace.tsv  | 2124  |
| DOMAINS/Chitin_synth_2.replace.tsv  | 2029  |

Check each domain, such as https://pfam.xfam.org/family/Epimerase. Some domains are not directly related to carbohydrate metabolism.

- ATP, AMP, GDP, CTP, and NADP associated
  - AAA_3
  - GDP_Man_Dehyd
  - CBS*
  - CTP_transf_like
- Binding to other small molecules
  - HAD* : enzymatic cleavage by nucleophilic substitution of carbon–halogen bonds (C–halogen) 卤素
  - F420: NADP oxidoreductase coenzyme F420-dependent
  - SKI: Shikimate kinase 莽草酸激酶 - ATP + shikimate <-> ADP + shikimate 3-phosphate
  - QPRTase_C: Quinolinate phosphoribosyl transferase, C-terminal domain 喹啉酸磷酸核糖转移酶
- Oxidoreductases
  - Epimerase: NADH dehydrogenase (ubiquinone) NADH脱氢酶泛醌 (辅酶Q)
  - ApbA: Ketopantoate reductase PanE/ApbA 酮戊酸还原酶 - redution alpha-ketopantoate -> D-(-)-pantoate (降解alpha-酮戊二酸为D-(-)-泛酸)
  - 3Beta_HSD: 3-beta-hydroxysteroid dehydrogenase 羟基化类固醇脱氢酶
  - NAD_Gly3P_dh_N: Glycerol-3-phosphate dehydrogenase 甘油-3-磷酸脱氢酶 - glycerol 3-phosphate -> dihydroxyacetone phosphate (3-磷酸甘油变为二羟丙酮磷酸)

We cannot `interproscan` all `200612` proteins at once.

```bash
cd /mnt/e/data/Pseudomonas

# All proteins appeared
cat pfam_domain.tsv | cut -f 2 | sort |
    tsv-filter --not-regex 1:'^AAA' |
    tsv-filter --not-regex 1:'^GDP' |
    tsv-filter --not-regex 1:'^CBS' |
    tsv-filter --not-regex 1:'^CTP' |
    tsv-filter --not-regex 1:'^HAD' |
    tsv-filter --not-regex 1:'^F420' |
    tsv-filter --not-regex 1:'^SKI' |
    tsv-filter --not-regex 1:'^QRPTase' |
    tsv-filter --not-regex 1:'^Epimerase' |
    tsv-filter --not-regex 1:'^ApbA' |
    tsv-filter --not-regex 1:'^3Beta_HSD' |
    tsv-filter --not-regex 1:'^NAD_Gly3P' \
    > domain.lst

wc -l < domain.lst
#  202

cat domain.lst |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        tsv-select -f 2 DOMAINS/{}.replace.tsv
    ' |
    sort -u \
    > DOMAINS/domains.tsv
# sort -u/--unique: with -c, check for strict ordering;
# without -c, output only the first of an equal run
# so this step actually sort and uniq them

# proteins have domains
wc -l < DOMAINS/domains.tsv
#  200612

# all proteins uniq among all strains
faops size PROTEINS/all.uniq.fa | wc -l
#  3944568

for domain in $(cat domain.lst); do
    echo 1>&2 "==> domain [${domain}]"

    tsv-join \
        DOMAINS/domains.tsv \
        --data-fields 1 \
        -f <(
            cat DOMAINS/${domain}.replace.tsv |
                perl -nla -e 'print qq{$F[1]\tO}'
        ) \
        --key-fields 1 \
        --append-fields 2 \
        --write-all "" \
        > DOMAINS/tmp.tsv

    mv DOMAINS/tmp.tsv DOMAINS/domains.tsv
done
# This step actually write all domains to a matrix
# Matrix raws represented domains, and cols meant proteins have domains
# "O" represented this protein existing that domain 

datamash check < DOMAINS/domains.tsv
#200612 lines, 203 fields, one more for proteins name fields

# Add header line
for domain in $(cat domain.lst); do
    echo "${domain}"
done |
    (echo -e "#name" && cat) |
    paste -s -d $'\t' - \
    > DOMAINS/header.tsv


cat DOMAINS/header.tsv DOMAINS/domains.tsv \
    > tmp.tsv && mv tmp.tsv DOMAINS/domains.tsv

tsv-join \
    PROTEINS/all.info.tsv \
    --data-fields 1 \
    -f DOMAINS/domains.tsv \
    --key-fields 1 \
    --append-fields 2-203 |
     keep-header -- sort -k1,1 \
    > tmp.tsv && mv tmp.tsv DOMAINS/domains.tsv

datamash check < DOMAINS/domains.tsv
#200613 lines, 206 fields
# another 3 fields added: stain, size and annotation

rm DOMAINS/header.tsv
```

### InterProScan

The previous `hmmsearch` steps were done to narrow down the number of target proteins so that as few proteins as possible would be passed to the following IPS steps.

InterProScan starts slowly, so dozens of proteins from one strain are submitted at once.

> - `interproscan.sh`
> 
> ```txt
> Options:
> -cpu/--cpu <CPU>: Optional, number of cores for interproscan
> -dp/--disable-precalc: Optional. Disables use of the precalculated match
>     lookup service.All match calculations will be run locally.
> -f/--formats <OUTPUT-FORMATS>: Optional, case-insensitive, comma separated
>     list of output formats. Supported formats are TSV, XML, JSON, GFF3, HTML and SVG.
> -i/--input <INPUT-FILE-PATH>: Optional, path to fasta file that should be
>     loaded on Master startup.
> -b/--output-file-base <OUTPUT-FILE-BASE>: Optional, base output filename (relative
>     or absolute path).
> ```
> 
> - `jq`
> 
> ```txt
> commandline JSON processor [version 1.6]
> 
> Usage:
> jq [options] <jq filter> [file...]
> jq [options] --args <jq filter> [strings...]
> jq [options] --jsonargs <jq filter> [JSON_TEXTS...]
> 
> jq is a tool for processing JSON inputs, aapplying the given filter to
> its JSON text inputs and producing the filter's results as JSON on
> standard output.
> 
> Options:
> -c: compact instead of pretty-printed output
> -r: output raw strings, not JSON texts
> ```

The goal is to use seq of proteins contained domains to search in each strain genome for similar ones.

```bash
cd /mnt/e/data/Pseudomonas

mkdir -p IPS

# extract wanted protein sequences in every strain
cat strains.lst |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        if [ $(({#} % 10)) -eq "0" ]; then
            >&2 printf "."
        fi

        mkdir -p IPS/{}

        cat PROTEINS/all.info.tsv |
            tsv-filter --str-eq 2:{} |
            cut -f 1 |
            grep -Fx -f <(cut -f 1 DOMAINS/domains.tsv | grep "^{}") \
            > IPS/{}/wanted.lst

        faops some PROTEINS/all.replace.fa IPS/{}/wanted.lst IPS/{}/{}.fa
    '
# {#} means sequence number of job to run
# So the meaning is every ten jobs finished will print a dot on screen
# each strain will give out a fasta file
# 1953 strains pass and be used for protein analysis

find IPS -maxdepth 1 -mindepth 1 -type d | wc -l
#  1952
# just equal to the number of strains

cat strains.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        cat IPS/{}/wanted.lst | wc -l
    ' |
    tsv-summarize --quantile 1:0.25,0.5,0.75
#  64      111     128
# totally 202 domains in domain.lst, and tsv-summarize could tell us quartile of domains

# scan proteins of each strain with InterProScan
# By default InterProScan uses 8 cpu cores 
mkdir -p split
split -a 4 -l 30 -d strains.lst split/
# split: output pieces of FILE to PREFIXaa, PREFIXab
# -a/--suffix-length=N: generate suffixes of length N (default 2).
# -d: use numeric suffixes starting at 0, not alphabetic 
# -a 4 -d: means the output file will have 4 suffix length, which means 0000, 0001...
# -l/--lines=NUMBER: put NUMBER lines/records per output file

# using protein sequences with domains to scan similar one
for f in $(find split -maxdepth 1 -type f -name "[0-9]*" | sort); do
    >&2 echo "==> IPS [${f}]"
    bsub -q mpi -n 24 -J "IPS-${f}" "
        cat ${f} |
            parallel --no-run-if-empty --linebuffer -k -j 6 '
                if [ -e IPS/{}/{}.tsv ]; then
                    >&2 echo {};
                    exit;
                fi

                interproscan.sh --cpu 4 -dp -f tsv,json -i IPS/{}/{}.fa --output-file-base IPS/{}/{}
            '
        "
done

rm -fr split output.*

find IPS -type f -name "*.json" | sort |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        if [ $(({#} % 10)) -eq "0" ]; then
            >&2 printf "."
        fi
        pigz -p 3 {}
    '
# -p/--processes n: allow up to n compression threads (default is the number of online processors) 

# IPS family
# Some proteins belong to more than one family. Only the best one is kept here
cat strains.lst |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        if [ $(({#} % 10)) -eq "0" ]; then
            >&2 printf "."
        fi
        gzip -dcf IPS/{}/{}.json.gz |
            jq .results |
            jq -r -c '\''
                .[] |
                .xref as $name |
                .matches[] |
                .signature.entry |
                select(.type == "FAMILY") |
                [$name[0].name, .accession, .description] |
                @tsv
            '\'' |
            tsv-uniq -f 1
    ' |
    (echo -e "#name\tfamily\tdescription" && cat) \
    > IPS/family.tsv

wc -l < IPS/family.tsv
#  135966
# proteins scanned by the method above

tsv-join \
    <(cut -f 1-4 DOMAINS/domains.tsv) \
    --data-fields 1 \
    -f IPS/family.tsv \
    --key-fields 1 \
    --append-fields 2-3 \
    --write-all "" |
    tsv-join \
        --data-fields 1 \
        -f DOMAINS/domains.tsv \
        --key-fields 1 \
        --append-fields 5-206 |
     keep-header -- sort -k1,1 \
    > IPS/predicts.tsv
```

InterPro provides functional analysis of proteins by classifying them into families and predicting domains and important sites. To classify proteins in this way, InterPro uses predictive models, known as signatures, provided by several different databases (referred to as member databases) that make up the InterPro consortium.

I currently know nothing about JSON, so the `jq` here could only be inferred from results.

`IPS/predicts.tsv` contains all proteins predicted by `InterProScan` with protein family name and description joined into `DOMAINS/domains.tsv`

## Search for gene families with uneven members

```bash
cd /mnt/e/data/Pseudomonas

cat IPS/predicts.tsv |
    tsv-filter -H --not-empty description |
    tsv-summarize -H -g strain,family,description --count |
    keep-header -- tsv-sort -k4,4n |
    tsv-summarize -H -g family,description --values 4 |
    keep-header -- tsv-sort -k2,2 |
    tsv-filter -H --str-in-fld 3:'|' |
    tr '|' '/' |
    keep-header -- perl -nla -F"\t" -e '
        my %seen;
        my @values = split "/", $F[2];
        next if scalar @values < 200; # The family is not widely distributed
        for my $v (@values) {
            $seen{$v}++;
        }
        next if keys %seen >= 5; # Too much variation in family members
        my @valid = grep { $seen{$_} > 50 } keys %seen;
        next if scalar @valid < 2; # At least 2 categories
        print join "\t", $F[0], $F[1], join "/", sort keys %seen;
    ' |
    tsv-filter -H --istr-not-in-fld 2:"probable" |
    tsv-filter -H --istr-not-in-fld 2:"putative" |
    tsv-filter -H --istr-not-in-fld 2:"Uncharacterised" \
    > variety.tsv
# find those proteins neither too common among all strains (too conservative)
# nor too variable in pro numbers among all strains (too diverse)

cat variety.tsv |
    mlr --itsv --omd cat
```

| family    | description                                               | count_values |
|-----------|-----------------------------------------------------------|--------------|
| IPR008183 | Aldose 1-/Glucose-6-phosphate 1-epimerase                 | 1/2/3/4      |
| IPR005841 | Alpha-D-phosphohexomutase superfamily                     | 1/2/3/4      |
| IPR014438 | Glucan biosynthesis protein MdoG/MdoD                     | 1/2/3/4      |
| IPR025532 | Glucose-6-phosphate 1-epimerase                           | 1/2          |
| IPR005999 | Glycerol kinase                                           | 1/2          |
| IPR011837 | Glycogen debranching enzyme, GlgX type                    | 1/2/3        |
| IPR000322 | Glycoside hydrolase family 31                             | 1/2/3/4      |
| IPR000811 | Glycosyl transferase, family 35                           | 1/2/3/4      |
| IPR035461 | GmhA/DiaA                                                 | 1/2          |
| IPR026040 | Hydroxypyruvate isomerase-like                            | 1/2/3/4      |
| IPR005501 | LamB/YcsF/PxpA-like                                       | 1/2/3/4      |
| IPR006415 | P-type ATPase, subfamily IIIB                             | 1/2/3        |
| IPR037950 | Peptidoglycan deacetylase PgdA-like                       | 1/2/3/4      |
| IPR004800 | Phosphosugar isomerase, KdsD/KpsF-type                    | 1/2/3        |
| IPR023853 | Poly-beta-1,6 N-acetyl-D-glucosamine synthase PgaC/IcaA   | 1/2          |
| IPR023854 | Poly-beta-1,6-N-acetyl-D-glucosamine N-deacetylase PgaB   | 1/2/3        |
| IPR004625 | Pyridoxine kinase                                         | 1/2/3        |
| IPR002347 | Short-chain dehydrogenase/reductase SDR                   | 1/2/3        |
| IPR017583 | Tagatose/fructose phosphokinase                           | 1/2/3        |
| IPR004730 | Transaldolase type 1                                      | 1/2/3/4      |
| IPR001585 | Transaldolase/Fructose-6-phosphate aldolase               | 1/2/3/4      |
| IPR005886 | UDP-glucose 4-epimerase                                   | 1/2/3        |
| IPR005593 | Xylulose 5-phosphate/Fructose 6-phosphate phosphoketolase | 1/2/3        |
| IPR005888 | dTDP-glucose 4,6-dehydratase                              | 1/2/3        |

### Within species

```bash
cd /mnt/e/data/Pseudomonas

for f in $(cat variety.tsv | tsv-select -f 1 | sed '1d' | sort); do
    cat IPS/predicts.tsv |
        tsv-select -f 1-3,5,6 |
        tsv-filter -H --str-eq family:"${f}" |
        tsv-join -d 2 \
            -f strains.taxon.tsv -k 1 \
            --append-fields 4 |
        tsv-summarize -g 6,2 --count |
        keep-header -- tsv-sort -k3,3n |
        tsv-summarize -g 1 --unique-values 3 |
        tsv-filter --str-in-fld 2:'|' |
        tsv-sort -k1,1 |
        tr '|' '/' \
        > tmp.tsv

    COUNT=$(cat tmp.tsv | wc -l)
    HAS_T=$(cat tmp.tsv | grep -E "aeruginosa|fluorescens|putida|syringae" | wc -l)
    HAS_ENOUGH=$(
        cat IPS/predicts.tsv |
            tsv-select -f 1-3,5,6 |
            tsv-filter -H --str-eq family:"${f}" |
            tsv-join -d 2 \
                -f strains.taxon.tsv -k 1 \
                --append-fields 4 |
            tsv-summarize -g 6,2 --count |
            tsv-filter --ge 3:2 |
            tsv-summarize -g 1 --count |
            tsv-filter --ge 2:5 |
            wc -l
    )
    if [[ ${COUNT} -ge "1" && ${COUNT} -le "10" && ${HAS_T} -ge "1" && ${HAS_ENOUGH} -ge "1" ]]; then
        cat tmp.tsv |
            sed "1 s/^/${f}\\t/" |
            sed "2,$ s/^/\\t/"
    fi

done |
    (echo -e "#family\tspecies\tcount" && cat) |
    mlr --itsv --omd cat
```

| #family   | species                    | count |
|-----------|----------------------------|-------|
| IPR000811 | Pseudomonas balearica      | 1/2   |
|           | Pseudomonas fluorescens    | 1/2   |
|           | Pseudomonas mandelii       | 1/2   |
|           | Pseudomonas stutzeri       | 1/2   |
|           | Pseudomonas synxantha      | 1/2   |
|           | Pseudomonas veronii        | 1/2   |
| IPR004730 | Alteromonas macleodii      | 1/2   |
|           | Alteromonas mediterranea   | 1/2   |
|           | Halomonas meridiana        | 1/2   |
|           | Halomonas titanicae        | 1/2/3 |
|           | Pseudomonas azotoformans   | 1/2   |
|           | Pseudomonas fluorescens    | 1/2   |
|           | Pseudomonas synxantha      | 1/2   |
| IPR004800 | Alteromonas macleodii      | 1/2   |
|           | Halomonas piezotolerans    | 1/2   |
|           | Pseudomonas atacamensis    | 1/2   |
|           | Pseudomonas citronellolis  | 1/2   |
|           | Pseudomonas fluorescens    | 1/2   |
|           | Pseudomonas fragi          | 1/2   |
|           | Pseudomonas putida         | 1/2   |
|           | Thalassolituus oleivorans  | 1/2   |
| IPR005593 | Pseudomonas aeruginosa     | 1/2   |
|           | Pseudomonas veronii        | 1/2   |
| IPR005999 | Alteromonas australica     | 1/2   |
|           | Pseudomonas aeruginosa     | 1/2   |
| IPR035461 | Alteromonas mediterranea   | 1/2   |
|           | Pseudomonas aeruginosa     | 1/2   |
|           | Pseudomonas lalkuanensis   | 1/2   |
| IPR037950 | Pseudomonas brassicacearum | 1/3   |
|           | Pseudomonas chlororaphis   | 1/2/3 |
|           | Pseudomonas fluorescens    | 1/2/3 |

### Among species

```bash
cd /mnt/e/data/Pseudomonas

for f in $(cat variety.tsv | tsv-select -f 1 | sed '1d' | sort); do
    cat IPS/predicts.tsv |
        tsv-select -f 1-3,5,6 |
        tsv-filter -H --str-eq family:"${f}" |
        tsv-join -d 2 \
            -f strains.taxon.tsv -k 1 \
            --append-fields 4 |
        tsv-summarize -g 6,2 --count |
        keep-header -- tsv-sort -k3,3n |
        tsv-summarize -g 1 --unique-values 3 |
        tsv-filter --str-not-in-fld 2:'|' |
        tsv-sort -k1,1 |
        tr '|' '/' \
        > tmp.tsv

    IS_1=$(cat tmp.tsv | tsv-filter --eq 2:1 | wc -l)
    IS_N=$(cat tmp.tsv | tsv-filter --ne 2:1 | wc -l)
    HAS_T=$(cat tmp.tsv | tsv-filter --ne 2:1 | grep -E "aeruginosa|chlororaphis|fluorescens|protegens|putida|syringae" | wc -l)

    if [[ ${IS_1} -ge "20" && ${IS_N} -ge "1" && ${IS_N} -le "10" && ${HAS_T} -ge "1" ]]; then
        printf "%s\t%d\t%d\t%d\n" $f $IS_1 $IS_N $HAS_T
        cat tmp.tsv | tsv-filter --ne 2:1 | grep -E "aeruginosa|chlororaphis|fluorescens|protegens|putida|syringae"
    fi
done
```

The goal is to find proteins with more than a copy within species (Pseudomonas), and better a copy among other species (not Pseudomonas).


### IPR005999 - Glycerol kinase

- IPR000577 - Carbohydrate kinase, FGGY
  - IPR005999 - Glycerol kinase
  - IPR006000 - Xylulokinase

```bash
cd /mnt/e/data/Pseudomonas

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR005999  |
    tsv-summarize -H -g annotation --count
#annotation      count
#glycerol kinase GlpK    1357
#glycerol kinase 7
#glycerol kinase (sn-glycerol-3-phosphate generating)    1

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR000577  |
    tsv-summarize -H -g annotation --count |
    tsv-filter -H --gt 2:5
#annotation      count
#hypothetical protein    9
#FGGY-family carbohydrate kinase 707
#glycerol kinase GlpK    839
#carbohydrate kinase     176
#xylulokinase    362
#L-fuculokinase  19
#autoinducer-2 kinase    21
#rhamnulokinase  15
#glycerol kinase 57
#pentose kinase  9
#sugar kinase    14
#ribulokinase    19
#FGGY family carbohydrate kinase 19

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR006000  |
    tsv-summarize -H -g annotation --count |
    tsv-filter -H --gt 2:5
#annotation      count
#xylulokinase    671

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR005999 |
    tsv-summarize -H -g size --count |
    keep-header -- tsv-sort -k1,1n |
    tsv-filter -H --ge count:10
#size    count
#493     10
#494     340
#495     22
#499     71
#500     32
#501     196
#502     255
#503     48
#504     38
#505     282
#506     11
#507     18

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR005999 |
    tsv-select -f 1-6 |
    tsv-join -d 2 \
        -f strains.taxon.tsv -k 1 \
        --append-fields 4 |
    tsv-filter --or --str-eq 7:"Pseudomonas aeruginosa" |
    tsv-summarize -g 7,2 --count |
    tsv-filter --gt 3:1 |
    tsv-summarize -g 1 --count
#Pseudomonas aeruginosa  195
```

All protein have the same structure FGGY_C and FGGY_N.

```bash
cd /mnt/e/data/Pseudomonas

mkdir -p /mnt/e/data/Pseudomonas/GlpK

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR005999 |
    datamash transpose |
    perl -nl -e '
        $row = $_;
        $row =~ s/\s//g;
        length($row) > 20 and print;
    ' |
    datamash transpose \
    > GlpK/GlpK.tsv

plotr tsv GlpK/GlpK.tsv --header

cat IPS/predicts.tsv |
    tsv-filter -H --or --str-eq family:IPR005999 --str-eq family:IPR000577 --str-eq family:IPR006000 |
    datamash transpose |
    perl -nl -e '
        $row = $_;
        $row =~ s/\s//g;
        length($row) > 20 and print;
    ' |
    datamash transpose \
    > GlpK/FGGY.tsv

cat DOMAINS/domains.tsv |
    tsv-filter -H --not-empty FGGY_C --not-empty FGGY_N |
    tsv-select -H -f 1-4,FGGY_C,FGGY_N \
    > GlpK/FGGY_N_C.tsv

wc -l GlpK/*.tsv
#  4320 GlpK/FGGY.tsv
#  5270 GlpK/FGGY_N_C.tsv
#  1364 GlpK/GlpK.tsv

for f in GlpK FGGY FGGY_N_C ; do
    >&2 echo "==> ${f}"

    faops some PROTEINS/all.replace.fa <(tsv-select -f 1 GlpK/${f}.tsv) stdout \
        > GlpK/${f}.fa

    muscle -in GlpK/${f}.fa -out GlpK/${f}.aln.fa

    FastTree GlpK/${f}.aln.fa > GlpK/${f}.aln.newick

    nw_reroot GlpK/${f}.aln.newick $(nw_labels GlpK/${f}.aln.newick | grep -E "Bac_subti|Sta_aure") |
        nw_order -c n - \
        > GlpK/${f}.reoot.newick

done
```

### IPR004800 - Phosphosugar isomerase, KdsD/KpsF-type

```bash
cd /mnt/e/data/Pseudomonas

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR004800 |
    tsv-summarize -H -g annotation --count
#annotation      count
#KpsF/GutQ family sugar-phosphate isomerase      1149
#D-arabinose 5-phosphate isomerase       8
#arabinose-5-phosphate isomerase 2
#carbohydrate isomerase  1
#D-arabinose 5-phosphate isomerase GutQ  1
#D-arabinose 5-phosphate isomerase KdsD  1
#arabinose-5-phosphate isomerase KdsD    389
#putative polysialic acid capsule expression protein     1
#hypothetical protein SF2731     1

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR004800 |
    tsv-summarize -H -g size --count |
    keep-header -- tsv-sort -k1,1n |
    tsv-filter -H --ge count:10
#size    count
#310     17
#315     13
#323     40
#324     452
#325     528
#326     422
#328     18
#339     16

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR004800 |
    tsv-select -f 1-6 |
    tsv-join -d 2 \
        -f strains.taxon.tsv -k 1 \
        --append-fields 4 |
    tsv-filter --or --str-eq 7:"Pseudomonas fluorescens" --str-eq 7:"Pseudomonas putida" |
    tsv-summarize -g 7,2 --count |
    tsv-filter --gt 3:1 |
    tsv-summarize -g 1 --count
#Pseudomonas fluorescens 2
#Pseudomonas putida      24
```

## InterProScan on all proteins of typical strains

Because the previous step is not able to find enough genes. It is obvious that another method should be taken as an alternative.

Using InterProScan to scan all proteins directly in typical strains of Pseudomonas species. Copy number of all proteins found will be counted and determined whether they are more than 1 copy in Pseudomonas aeruginosa.

```bash
cd /mnt/e/data/Pseudomonas

# the number of all proteins in PAO1
faops size ASSEMBLY/Pseudom_aeru_PAO1/*_protein.faa.gz |
    wc -l
#  5572

# the length of all proteins added
faops size ASSEMBLY/Pseudom_aeru_PAO1/*_protein.faa.gz |
    tsv-summarize --sum 2
#1858983

mkdir -p STRAINS

# 12 RefSeq strains of all Pseudomonas
for S in \
    Pseudom_aeru_PAO1 \
    Pseudom_puti_KT2440_GCF_000007565_2 \
    Pseudom_chl_aureofaciens_30_84_GCF_000281915_1 \
    Pseudom_entomophi_L48_GCF_000026105_1 \
    Pseudom_fluo_SBW25_GCF_000009225_2 \
    Pseudom_prot_Pf_5_GCF_000012265_1 \
    Pseudom_sav_pv_phaseolicola_1448A_GCF_000012205_1 \
    Pseudom_stu_A1501_GCF_000013785_1 \
    Pseudom_syr_pv_syringae_B728a_GCF_000012245_1 \
    Pseudom_aeru_UCBPP_PA14_GCF_000014625_1 \
    Pseudom_aeru_PA7_GCF_000017205_1 \
    Pseudom_aeru_LESB58_GCF_000026645_1 \
    ; do
    echo ${S}
done \
    > typical.lst

for S in $(cat typical.lst); do
    mkdir -p STRAINS/${S}
    faops split-about ASSEMBLY/${S}/*_protein.faa.gz 200000 STRAINS/${S}/
done
# faops split-about: split an fa file into several files of about approx_size bytes each by record
# faops split-about [options] <in.fa> <approx_size> <outdir>
# it was split by file size instead of number of lines

# using interproscan to search proteins in every RefSeq strains
# interproscan would consume a lot of computing power, so files were splited and executed on HPCC
for S in $(cat typical.lst); do
    for f in $(find STRAINS/${S}/ -maxdepth 1 -type f -name "[0-9]*.fa" | sort); do
        >&2 echo "==> ${f}"
        if [ -e ${f}.tsv ]; then
            >&2 echo ${f}
            continue
        fi

        bsub -q mpi -n 24 -J "${f}" "
            interproscan.sh --cpu 24 -dp -f tsv,json -i ${f} --output-file-base ${f}
        "
    done
done

find STRAINS -type f -name "*.json" | sort |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        if [ $(({#} % 10)) -eq "0" ]; then
            >&2 printf "."
        fi
        pigz -p 3 {}
    '

# same protein may have multiple families
for S in $(cat typical.lst); do
    for f in $(find STRAINS/${S} -maxdepth 1 -type f -name "[0-9]*.json.gz" | sort); do
        >&2 echo "==> ${f}"
        gzip -dcf ${f} |
            jq .results |
            jq -r -c '
                .[] |
                .xref as $name |
                .matches[] |
                .signature.entry |
                select(.type == "FAMILY") |
                [$name[0].name, .accession, .description] |
                @tsv
            ' |
            tsv-uniq
    done \
        > STRAINS/${S}/family.tsv
done

COUNT=
for S in $(cat typical.lst); do
    if [ ! -s STRAINS/${S}/family.tsv ]; then
        continue
    fi
    cat STRAINS/${S}/family.tsv |
        tsv-summarize -g 2,3 --count \
        > STRAINS/${S}/family-count.tsv

    COUNT=$((COUNT + 1))
done
echo $COUNT
#12
# $COUNT is the variable used for loop control
# it is also a proxy of all strains with protein families detected
# so the $COUNT equals to 12 also represents all strains with pro_families
# it equals to the pseudomonas RefSeq strains number

# families in all strains
for S in $(cat typical.lst); do
    cat STRAINS/${S}/family-count.tsv
done |
    tsv-summarize -g 1,2 --count |
    tsv-filter -H --istr-not-in-fld 2:"probable" |
    tsv-filter -H --istr-not-in-fld 2:"putative" |
    tsv-filter -H --istr-not-in-fld 2:"Uncharacterised" |
    tsv-filter -H --istr-not-in-fld 2:" DUF" |
    tsv-filter --ge 3:$COUNT \
    > STRAINS/universal.tsv
# for loop will cat all family-count into one
# tsv-summarize will count according to the family name - so it represents directly to copy number of 12 strains
# using tsv-filter to exclude those unknown proteins, meanwhile keep those proteins more than 12 strains existed
# --ge 3:$COUNT: means to get those genes with 1 or more copies in each strain
# this step actually ensure that every strain has at least a copy

# All other strains should have only 1 family member
cp STRAINS/universal.tsv STRAINS/family-1.tsv
for S in $(cat typical.lst | grep -v "_aeru_"); do
    if [ ! -s STRAINS/${S}/family-count.tsv ]; then
        continue
    fi
    cat STRAINS/${S}/family-count.tsv |
        tsv-join -k 1 -f STRAINS/family-1.tsv |
        tsv-filter --eq 3:1 \
        > STRAINS/family-tmp.tsv

    mv STRAINS/family-tmp.tsv STRAINS/family-1.tsv
done
# tsv-join to exclude those strains not in universal.tsv, which means copy number at least greater than 12
# so this step will find those only 1 family (proteins) in Pseudomonas strains (not aeruginosa)
# then tsv-filter will keep those only have 1 copy

# All P_aeru strains should have multiple family members
cp STRAINS/family-1.tsv STRAINS/family-n.tsv
for S in $(cat typical.lst | grep "_aeru_"); do
    if [ ! -s STRAINS/${S}/family-count.tsv ]; then
        continue
    fi
    cat STRAINS/${S}/family-count.tsv |
        tsv-join -k 1 -f STRAINS/family-n.tsv |
        tsv-filter --gt 3:1 \
        > STRAINS/family-tmp.tsv

    wc -l < STRAINS/family-tmp.tsv
    mv STRAINS/family-tmp.tsv STRAINS/family-n.tsv
done
#  18
#  16
#  14
#  14
# four numbers means the family-n.tsv was filtered four times as the P_aeru has four strains
# after that, 14 target genes were left for the next step

wc -l STRAINS/Pseudom_aeru_PAO1/family.tsv STRAINS/universal.tsv STRAINS/family-1.tsv STRAINS/family-n.tsv
#  4084 STRAINS/Pseudom_aeru_PAO1/family.tsv
#  1567 STRAINS/universal.tsv
#   972 STRAINS/family-1.tsv
#    14 STRAINS/family-n.tsv

cat STRAINS/family-n.tsv |
    tsv-select -f 1,2 |cc
    (echo -e "#family\tcount" && cat) |
    mlr --itsv --omd cat
```

| #family   | count                                                         |
|-----------|---------------------------------------------------------------|
| IPR014311 | Guanine deaminase                                             |
| IPR001404 | Heat shock protein Hsp90 family                               |
| IPR005999 | Glycerol kinase                                               |
| IPR000813 | 7Fe ferredoxin                                                |
| IPR011757 | Lytic transglycosylase MltB                                   |
| IPR007416 | YggL 50S ribosome-binding protein                             |
| IPR004361 | Glyoxalase I                                                  |
| IPR024922 | Rubredoxin                                                    |
| IPR001353 | Proteasome, subunit alpha/beta                                |
| IPR002307 | Tyrosine-tRNA ligase                                          |
| IPR024088 | Tyrosine-tRNA ligase, bacterial-type                          |
| IPR037532 | Peptidoglycan D,D-transpeptidase FtsI                         |
| IPR003672 | CobN/magnesium chelatase                                      |
| IPR004685 | Branched-chain amino acid transport system II carrier protein |

### IPR007416 - YggL 50S ribosome-binding protein

- Pfam: PF04320 - YggL_50S_bp
- PANTHER: PTHR38778 - CYTOPLASMIC PROTEIN-RELATED (PTHR38778)

```bash
cd /mnt/e/data/Pseudomonas

cat STRAINS/Pseudom_aeru_PAO1/*.tsv |
    grep "IPR007416"
# this step will show you how to get the above two accessions

mkdir -p YggL/HMM

curl -L https://pfam.xfam.org/family/PF04320/hmm > YggL/HMM/YggL_50S_bp.hmm
curl -L www.pantherdb.org/panther/exportHmm.jsp?acc=PTHR38778 > YggL/HMM/PTHR38778.hmm

# Ribosomal protein L10 and S8
curl -L https://pfam.xfam.org/family/PF00466/hmm > YggL/HMM/Ribosomal_L10.hmm
curl -L https://pfam.xfam.org/family/PF00410/hmm > YggL/HMM/Ribosomal_S8.hmm

E_VALUE=1e-20
for domain in YggL_50S_bp PTHR38778 Ribosomal_L10 Ribosomal_S8 ; do
    >&2 echo "==> domain [${domain}]"

    if [ -e YggL/${domain}.replace.tsv ]; then
        continue;
    fi

    for ORDER in $(cat order.lst); do
        >&2 echo "==> ORDER [${ORDER}]"

        cat taxon/${ORDER} |
            parallel --no-run-if-empty --linebuffer -k -j 8 "
                gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                    hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw YggL/HMM/${domain}.hmm - |
                    grep '>>' |
                    perl -nl -e '
                        m{>>\s+(\S+)} or next;
                        \$n = \$1;
                        \$s = \$n;
                        \$s =~ s/\.\d+//;
                        printf qq{%s\t%s_%s\n}, \$n, {}, \$s;
                    '
            "
    done \
        > YggL/${domain}.replace.tsv

    >&2 echo
done
# basic operations on hmmsearch, go check the above
# this E_VALUE will help us find domains most likely as usual

tsv-join YggL/YggL_50S_bp.replace.tsv \
    -f YggL/PTHR38778.replace.tsv \
    > YggL/YggL.replace.tsv
# combine two files included domains from different database

wc -l YggL/*.tsv
#  1559 YggL/PTHR38778.replace.tsv
#  1948 YggL/Ribosomal_L10.replace.tsv
#  1950 YggL/Ribosomal_S8.replace.tsv
#  1559 YggL/YggL.replace.tsv
#  1564 YggL/YggL_50S_bp.replace.tsv

# also extract all pro_seq
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 YggL/YggL.replace.tsv) YggL/YggL.fa

# muscle for multi_seq alignment
muscle -in YggL/YggL.fa -out YggL/YggL.aln.fa

FastTree YggL/YggL.aln.fa > YggL/YggL.aln.newick

nw_reroot YggL/YggL.aln.newick $(nw_labels YggL/YggL.aln.newick | grep -E "Bac_subti|Sta_aure") |
    nw_order -c n - \
    > YggL/YggL.reoot.newick
```

## Collect CDS

Getting all cds, maybe needed later.

### `all.cds.fa`

```bash
cd /mnt/e/data/Pseduomonas

mkdir -p CDS

find ASSEMBLY -type f -name "*_cds_from_genomic.fna.gz" |
    wc -l
#  1958

# sed script converting from Contigs to Strain
for ORDER in $(cat order.lst); do
    echo 1>&2 "==> ORDER [${ORDER}]"

    for STRAIN in $(cat taxon/${ORDER}); do
        find ASSEMBLY/${STRAIN} -type f -name "*_genomic.fna.gz" |
            grep -v "_from_" |
            xargs gzip -dcf |
            grep '^>' |
            cut -d' ' -f 1 |
            sed 's/>//' |
            xargs -I{} echo -e "{}\t${STRAIN}"
    done
done \
    > CDS/contigs_to_strain.tsv

cat CDS/contigs_to_strain.tsv |
    perl -nla -e '
        print q{s/^>} . quotemeta($F[0]) . q{/>} . quotemeta($F[1]) . q{/g;};
    ' \
    > CDS/sed.script

wc -l < CDS/sed.script
# 3893

for ORDER in $(cat order.lst); do
    echo 1>&2 "==> ORDER [${ORDER}]"

    for STRAIN in $(cat taxon/${ORDER}); do
        gzip -dcf ASSEMBLY/${STRAIN}/*_cds_from_genomic.fna.gz
    done
done |
    perl -nl -e 's/^>lcl\|/>/g; print' |
    perl -nl -e 's/\s+\[.+?\]//g; print' \
    > CDS/all.cds.fa
```

### `YggL.cds.fa`

```bash
cd /mnt/e/data/Pseudomonas

for domain in YggL Ribosomal_L10 Ribosomal_S8 ; do
    cat CDS/all.cds.fa |
        grep '>' |
        grep -F -f <( cat YggL/${domain}.replace.tsv | cut -f 1 ) |
        sed 's/^>//' \
        > CDS/${domain}.lst
done

for domain in YggL Ribosomal_L10 Ribosomal_S8 ; do
    faops order CDS/all.cds.fa CDS/${domain}.lst stdout |
        sed -f CDS/sed.script \
        > CDS/${domain}.cds.fa
done
```
