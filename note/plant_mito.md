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
