# Plant mitochondrion genomes

This script was originated from [plant_mito.md](https://github.com/wang-q/withncbi/blob/master/taxon/plant_mito.md). I noted and recorded my study.

## Preparation

### Download software

```bash
brew install miller librsvg
brew install mash newick_utils
brew install wang-q/tap/nwr wang-q/tap/tsv-utils
```

Install [egaz](https://github.com/wang-q/App-Egaz#installation)

### Update taxonomy database

```bash
rm -rf ~/.nwr
nwr download
nwr txdb
```

## Scrap id and acc from NCBI

Use `taxon/gb_taxon_locus.pl` to extract info from RefSeq files.

```bash
mkdir -p /mnt/e/data/mito/GENOMES
cd /mnt/e/data/mito/GENOMES

wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.1.genomic.gbff.gz
wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.2.genomic.gbff.gz

gzip -dcf mitochondrion.*.genomic.gbff.gz > genomic.gbff

perl ~/Scripts/withncbi/taxon/gb_taxon_locus.pl genomic.gbff > refseq_id_seq.csv
#There are [13200] sequences.
#There are [13200] valid sequences.

rm genomic.gbff
#rm *.genomic.gbff.gz

cat refseq_id_seq.csv | grep -v "^#" | wc -l
#13200

# combine
cat refseq_id_seq.csv |
    sort -u | # duplicated id-seq pair
    sort -t, -k1,1 |
    mlr --icsv --otsv cat \
    > id_seq.tsv
# sort
# -u/--unique: with -c, check for strict ordering; without -c, output only the first of an equal run
# -t/--field-separator=SEP: use SEP instead of non-blank to blank transition

cat id_seq.tsv | grep -v "^#" | wc -l
#13200
```

### Restrict taxonomy ids to green plants

```txt
Eukaryota (2759)
    Viridiplantae (33090) # 绿色植物界
        Chlorophyta (3041) # 绿藻门
        Streptophyta (35493) # 链型植物门
```

All living green plants belong to the major phylums including Streptophyta and Chlorophyta.

```bash
cd /mnt/e/data/mito/GENOMES

# Viridiplantae 33090
echo -e '#tax_id\taccession' > plant_id_seq.tsv # add headline
cat id_seq.tsv |
    nwr restrict 33090 -f stdin -c 1 \
    >> plant_id_seq.tsv

cat plant_id_seq.tsv | grep -v "^#" | wc -l
#414

# find repeated tax_id
cat plant_id_seq.tsv |
    cut -f 1 |
    sort -n |
    tsv-uniq --number --repeated |
    nwr append stdin
#3659    2       Cucumis sativus
#3659    3       Cucumis sativus
#3708    2       Brassica napus
#51329   2       Polytomella parva
#351366  2       Polytomella piriformis
# tsv-uniq --repeated: output only lines that are repeated (based on the key)
```

Look inside `plant_id_seq.tsv` and remove redundancies

```bash
# Cucumis sativus has 3 chromosomes
# Brassica napus linear plasmid NC_004946
# Polytomella parva has 2 chromosomes
# Polytomella piriformis (includes Polytomella sp. Pringsheim 63-10) has 2 chromosomes

sed -i".bak" "/NC_004946$/d" plant_id_seq.tsv #  napus
# .bak: often means filename extension commonly useBrassicad to signify a backup copy of a file
```

## Add lineage information

Give ids better shapes for manually checking and automatic filtering.

If you sure, you can add or delete lines and contents in `CHECKME.tsv`.

```bash
mkdir -p /mnt/e/data/mito/summary
cd /mnt/e/data/mito/summary

# generate a TSV file for manually checking
cat ../GENOMES/plant_id_seq.tsv |
    nwr append stdin |
    nwr append stdin -r species -r genus -r family -r order -r class -r phylum |
    keep-header -- sort -k9,9 -k8,8 -k7,7 -k6,6 -k5,5 \
    > CHECKME.tsv
```
