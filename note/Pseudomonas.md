# *Pseudomonas* HGT

This markdown is to record my study of teacher's [Pseudomonas.md](https://github.com/wang-q/withncbi/blob/master/pop/Pseudomonas.md).

Before getting start, the Taxonomy is very important for extracting those strains from all ranks.

The basic taxonomy ranks are here:

| Domain | Kingdom | Phylum | Class | Order | Family | Genus | Species |
|--------|---------|--------|-------|-------|--------|-------|---------|
| 域      | 界       | 门      | 纲     | 目     | 科      | 属     | 种       |

## Software

- PPanGGOLiN

Installation steps can be found [here](https://github.com/wang-q/dotfiles/blob/master/others.sh) 

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
```

RefSeq:

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

## Download all assemblies

- Get all assemblies info

```bash
cd /mnt/e/data/Pseudomonas

# Pseudomonas aeruginosa PAO1 is in the reference list
cat reference.tsv |
    tsv-select -H -f organism_name,species,genus,ftp_path,assembly_level \
    > raw.tsv
cat raw.tsv | wc -l
#16

# Species with 2 or more genomes were retained
SPECIES=$(
    cat species.count.tsv |
        tsv-filter -H --ge CHR:2 |
        sed '1d' |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

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

cat raw.tsv | wc -l
#1666

# The SQLite || operator allows you to concatenate 2 or more strings together.
# This step finally got all strains' assemblies from species.count.tsv

# Also includes representative strains of Gammaproteobacteria.
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
#2095

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
# A double dash (--) delimits the command, similar to how the pipe
# operator (|) delimits commands. Examples:
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

find biosample -name "SAM*.txt" | wc -l
# 1957

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
#1958 ASSEMBLY/n50.pass.csv
#1959 ASSEMBLY/n50.tsv
#3917 total

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
#1959 ASSEMBLY/Pseudomonas.assembly.collect.csv
#1953 ASSEMBLY/Pseudomonas.assembly.pass.csv
#3912 total
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

# Group by order
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
# All OK

rm *.tmp
```

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
#1952

cat dist.tsv | wc -l
#1904176
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
#3808352
# actually swap col1 and col2 and then using cat to combine them

cut -f 1 dist.tsv | tsv-uniq | parallel -j 1 --keep-order 'echo -e "{}\t{}\t0"' | wc -l
# other 1951 strains dist to themselves - 0

cat dist_full.tsv | wc -l
#3810303
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
# 391

ppanggolin workflow --anno pangenome/Pseudom_aeru_1.gbff.list --cpu 8 -o pangenome/Pseudom_aeru
# ppanggolin workflow
# .gbff.list example could be seen above
```

## Collect proteins

- `all.pro.fa`

Merge all proteins to a file and using accession to get uniq one.

```bash
cd /mnt/e/data/Pseudomonas

mkdir -p PROTEINS

find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
    sort |
    grep 'ASSEMBLY/' |
    wc -l
#1958

find ASSEMBLY -type f -name "*_protein.faa.gz" |
    wc -l
#1958

cat strains.lst |
    wc -l
#1952

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
#8754303

cat PROTEINS/all.pro.fa |
    grep "^>" |
    tsv-uniq |
    wc -l
#3995205

# annotations may be different
cat PROTEINS/all.uniq.fa |
    grep "^>" |
    wc -l
#3944568
# which means the same accession annotated differently 

# ribonuclease
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
#8754303

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
#8754303

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
#8754304
```

## Phylogenetics with 40 single-copy genes

### Find correspoding proteins by `hmmsearch`

- Download HMM models as described in `hmm/README.md`
- The `E_VALUE` was manually adjusted to 1e-20 to reach a balance between sensitivity and speciality.

```bash
E_VALUE=1e-20

cd /mnt/e/data/Pseudomonas

## example
#gzip -dcf ASSEMBLY/Ac_axa_ATCC_25176/*_protein.faa.gz |
#    hmmsearch -E 1e-20 --domE 1e-20 --noali --notextw ~/data/HMM/scg40/bacteria_and_archaea_dir/BA00001.hmm - |
#    grep '>>' |
#    perl -nl -e '/>>\s+(\S+)/ and print $1'

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
```
