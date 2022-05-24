# *Pseudomonas* HGT

This markdown is to record my study of teacher's [Pseudomonas.md](https://github.com/wang-q/withncbi/blob/master/pop/Pseudomonas.md).

## Strain info

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
| ---------------- | ----- |
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
| ---------------- | ----- |
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
| ------- | ---------------------------------------- | ---------------- |
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
# \b allows you to perform a “whole words only” search using a regular expression in the form of \bword\b

# nwr:
# nwr member:
# output format: #tax_id sci_name        rank    division
# nwr append:
# append fields of higher ranks to a TSV file
# * If `--rank` is empty, the scientific name will be appended.

# tsv-sort:
# tsv-sort runs 'sort' using TAB as the field delimiter, all arguments are forwared to 'sort'
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
| ---------- | ------------------------------- | ---- | --- |
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

# other bacteria genura
GENUS=$(
    nwr member Bacteria -r genus |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        sed '1d' |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//' # remove the last ','
)

# extract RefSeq from all genera
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
#16 - only 16 reference genome from all genera

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
# if 
```

RefSeq:

| #tax_id | organism_name                                                    | phylum                               |
| ------- | ---------------------------------------------------------------- | ------------------------------------ |
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
