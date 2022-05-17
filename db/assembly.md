# NCBI Assembly Reports

Download date: 2022-5-16

This script was originated from my teacher [wangq](https://github.com/wang-q/nwr/blob/master/doc/assembly.md).

## Preparations

Software you need.

```bash
brew install wangq/tap/nwr
brew install wangq/tap/tsv-utils
brew install sqlite
brew install miller

# Download `taxdump.tar.gz` and assembly reports
nwr download

# Init the taxonomy database
nwr txdb

# Init the assembly database
nwr ardb
nwr ardb --genbank
```

## NCBI ASSEMBLY

- assembly_level

```bash
for C in refseq genbank; do
    cat ~/.nwr/assembly_summary_${C}.txt |
        sed '1d' |
        tsv-summarize -H -g assembly_level,genome_rep --count |
        keep-header -- sort |
        mlr --itsv --omd cat
    
    echo -e "\nTable: ${C}\n\n"
done
```

| assembly_level  | genome_rep | count  |
| --------------- | ---------- | ------ |
| Chromosome      | Full       | 4812   |
| Chromosome      | Partial    | 433    |
| Complete Genome | Full       | 38618  |
| Complete Genome | Partial    | 2      |
| Contig          | Full       | 133752 |
| Contig          | Partial    | 1      |
| Scaffold        | Full       | 83752  |
| Scaffold        | Partial    | 28     |

Table: refseq

| assembly_level  | genome_rep | count   |
| --------------- | ---------- | ------- |
| Chromosome      | Full       | 9234    |
| Chromosome      | Partial    | 2045    |
| Complete Genome | Full       | 78166   |
| Complete Genome | Partial    | 146     |
| Contig          | Full       | 1004474 |
| Contig          | Partial    | 845     |
| Scaffold        | Full       | 177874  |
| Scaffold        | Partial    | 367     |

Table: genbank

## Count qualified assemblies of Prokaryote groups

```bash
mkdir -p /mnt/d/data/bacteria/ASSEMLY
cd /mnt/d/data/bacteria

echo -e "GROUP_NAME\tComplete Genome\tChromosome\tScaffold\tContig" \
    > groups.tsv

for item in Bacteria Archaea ; do
    PHYLUM=$(
        nwr member ${item} -r phylum |
            grep -v -i "Candidatus " |
            grep -v -i "candidate " |
            sed '1d' |
            cut -f 2 |
            sort
    )
    echo -e "$item\t\t\t\t"

    for P in $PHYLUM; do
        GENUS=$(
            nwr member ${P} -r genus |
                grep -v -i "Candidatus " |
                grep -v -i "candidate " |
                sed '1d' |
                cut -f 1 |
                tr "\n" "," |
                sed 's/,$/\)/' |
                sed 's/^/\(/'
        )
        
        if [[ ${#GENUS} -lt 3 ]]; then
            >&2 echo $P has no genera
            continue 
        fi
    
        printf "$P\t"
    
        for L in 'Complete Genome' 'Chromosome' 'Scaffold' 'Contig'; do
            echo "
                SELECT
                    COUNT(*)
                FROM ar
                WHERE 1=1
                    AND genus_id IN $GENUS
                    AND assembly_level IN ('$L')
                " |
                sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
        done |
        tr "\n" "\t" |
        sed 's/\t$//'
        
        echo;
    done
done  \
    >> groups.tsv

cat groups.tsv |
    mlr --itsv --omd cat
```

| GROUP_NAME            | Complete Genome | Chromosome | Scaffold | Contig |
| --------------------- | --------------- | ---------- | -------- | ------ |
| Bacteria              |                 |            |          |        |
| Abditibacteriota      | 0               | 0          | 0        | 1      |
| Acidobacteria         | 21              | 1          | 27       | 24     |
| Actinobacteria        | 2625            | 526        | 9421     | 10466  |
| Aquificae             | 14              | 2          | 8        | 9      |
| Armatimonadetes       | 1               | 1          | 4        | 1      |
| Atribacterota         | 1               | 0          | 0        | 0      |
| Bacteroidetes         | 968             | 176        | 2738     | 3573   |
| Balneolaeota          | 0               | 0          | 4        | 16     |
| Caldiserica           | 1               | 0          | 0        | 0      |
| Calditrichaeota       | 1               | 1          | 0        | 0      |
| Chlamydiae            | 189             | 83         | 52       | 102    |
| Chlorobi              | 13              | 0          | 6        | 11     |
| Chloroflexi           | 4               | 0          | 4        | 4      |
| Chrysiogenetes        | 2               | 0          | 2        | 0      |
| Coprothermobacterota  | 1               | 0          | 1        | 2      |
| Cyanobacteria         | 191             | 45         | 257      | 431    |
| Deferribacteres       | 5               | 0          | 3        | 7      |
| Deinococcus-Thermus   | 71              | 3          | 58       | 134    |
| Dictyoglomi           | 2               | 0          | 0        | 1      |
| Elusimicrobia         | 2               | 0          | 0        | 1      |
| Fibrobacteres         | 2               | 0          | 10       | 28     |
| Firmicutes            | 5954            | 897        | 27671    | 35457  |
| Fusobacteria          | 78              | 5          | 102      | 135    |
| Gemmatimonadetes      | 4               | 0          | 2        | 1      |
| Ignavibacteriae       | 2               | 0          | 0        | 0      |
| Kiritimatiellaeota    | 2               | 0          | 0        | 2      |
| Lentisphaerae         | 0               | 0          | 2        | 4      |
| Nitrospinae           | 0               | 0          | 1        | 2      |
| Nitrospirae           | 9               | 0          | 3        | 10     |
| Planctomycetes        | 54              | 26         | 37       | 49     |
| Proteobacteria        | 15377           | 2189       | 41410    | 80746  |
| Rhodothermaeota       | 0               | 0          | 3        | 3      |
| Spirochaetes          | 199             | 139        | 267      | 833    |
| Synergistetes         | 6               | 4          | 10       | 20     |
| Tenericutes           | 426             | 18         | 161      | 406    |
| Thermodesulfobacteria | 7               | 0          | 4        | 5      |
| Thermotogae           | 41              | 1          | 32       | 38     |
| Verrucomicrobia       | 112             | 7          | 157      | 92     |
| Archaea               |                 |            |          |        |
| Crenarchaeota         | 93              | 9          | 10       | 77     |
| Euryarchaeota         | 292             | 9          | 248      | 411    |
| Nanoarchaeota         | 0               | 0          | 0        | 0      |
| Thaumarchaeota        | 10              | 0          | 4        | 4      |

Table: refseq - Prokaryotes

## Update old species genomes

Extract species information from old [ASSEMBLY](https://github.com/wang-q/withncbi/tree/master/db) steps.

The goal of this step is to get all strains good assembly from the species information.

### Strain info

```bash
cd /mnt/d/data/bacteria/ASSEMBLY

cat target.assembly.collect.csv | 
    tsv-select -H -f Organism_name |
    sed '1d' |
    tsv-select -d " " -f 1,2
    > ../organism_name.tsv

cd ..

# all species info in nwr
cat organism_name.tsv |
    parallel --line-buffer '
        nwr member {} |
        grep -v '^\#' >> organism_name.all.tsv
    '
sed -i '1i#tax_id\tsci_name\trank\tdivision' organism_name.all.tsv

cat organism_name.all.tsv | 
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat
```

| rank               | count |
| ------------------ | ----- |
| species            | 186   |
| strain             | 24562 |
| subspecies         | 166   |
| no rank            | 1900  |
| biotype            | 7     |
| isolate            | 124   |
| serogroup          | 134   |
| serotype           | 234   |
| forma specialis    | 437   |
| Salmonella bongori | 1     |
| Bacteria           | 1     |
| varietas           | 1     |

```bash
cat organism_name.tsv |
    parallel --line-buffer '
        nwr member {} -r species -r strain |
        grep -v '^\#' >> organism_name.strain.tsv
    '
sed -i '1i#tax_id\tsci_name\trank\tdivision' organism_name.strain.tsv

cat organism_name.strain.tsv | 
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat
```

| rank                  | count |
| --------------------- | ----- |
| species               | 186   |
| strain                | 24562 |
| Gluconobacter oxydans | 1     |
| Bacteria              | 1     |

### Species with assemblies

Extract all those strains from 186 species.

```bash
cd /mnt/d/data/bacteria

SPECIES=$(
    cat organism_name.strain.tsv |
        tsv-filter -H --str-eq rank:species |
        tsv-select -f 1 |
        sed '1d' |
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


```
