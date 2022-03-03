#  Get data from NCBI

Using rsync or aspera

```bash
rsync -avP ftp.ncbi.nlm.nih.gov::genomes/GENOME_REPORTS/ \
    /mnt/e/data/NCBI/genomes/GENOME_REPORTS/
```

```bash
rsync -avP ftp.ncbi.nlm.nih.gov::genomes/ASSEMBLY_REPORTS/ \
    --exclude=".tmp" \
    --exclude=".old" \
    /mnt/e/data/NCBI/genomes/ASSEMBLY_REPORTS/
```

```bash
rsync -avP ftp.ncbi.nlm.nih.gov::bioproject/ \
    --exclude="*.xml" \
    /mnt/e/data/NCBI/bioproject/
```

```bash
rsync -avP ftp.ncbi.nlm.nih.gov::pub/taxonomy/ \
    --exclude=".tmp" \
    --exclude=".old" \
    --exclude="*.Z" \
    --exclude="taxdump_archive" \
    --exclude="new_taxdump" \
    --exclude="accession2taxid" \
    --exclude="gi_taxid_*" \
    /mnt/e/data/NCBI/taxonomy/

rm -fr /mnt/e/data/NCBI/taxdmp
mkdir -p /mnt/e/data/NCBI/taxdmp
tar xvfz /mnt/e/data/NCBI/taxonomy/taxdump.tar.gz -C /mnt/e/data/NCBI/taxdmp
```

Because my WSL2 broke for GFW, I just can't install brew for manipulating apps.

#  Download Staphylococcus genomes

##  Genomes from bacteria_gr

### Choose strains from genus Staphylococcus

```bash
mkdir -p ~/data/bacteria/summary
cd ~/data/bacteria/summary
mysql -ualignDB -palignDB gr_prok -e '
	SELECT taxonomy_id `#strain_taxonomy_id`,
           organism_name `strain`,
           species `species`,
           genus `genus`,
           subgroup `subgroup`,
           `code`,
           chr
    FROM gr
	WHERE 1=1
	AND genus LIKE "%Staphylococcus%"
    AND (status LIKE "%Complete%"
    	OR status LIKE "%Chromosome%")
    AND species NOT LIKE "%Candidatus%"
    AND taxonomy_id != species_id
    AND organism_name NOT LIKE "%,%"
    AND (chr IS NOT NULL OR LENGTH(CHR) > 0)
    AND genus IS NOT NULL
    AND (LENGTH(wgs) = 0 OR wgs IS NULL)
	' | 
	grep -v "^#" > GENUS_Sta.tsv

cat GENUS_Sta.tsv |
	mlr --itsv --ocsv cat > GENUS_Sta.csv
```

### Init genome report database

```bash
echo '#strain_taxonomy_id,strain,species,genus,subgroup,code,accession,abbr' > ABBR.csv.tmp
cat GENUS_Sta.csv | grep -v '^#' |
	perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "2,3,4" -s "," -m 0 --shortsub |
	sort -t',' -k5,5 -k4,4 -k3,3 -k6,6 >> ABBR.csv.tmp
	
cat ABBR.csv.tmp |
	grep -v "^#" |
    cut -d, -f 3 |
    uniq -c |
    sort -nr |
    head -n 10

cat ABBR.csv.tmp |
	perl -nla -F"," -e '
		if (
			$F[4] eq q{Staphylococcus}
			) {
				$F[4] =~ /^NZ_/ and next;
			}
			
			print;
	' > ABBR.csv

cat ABBR.csv |
	grep -v "^#" |
    cut -d, -f 3 |
    uniq -c |
    sort -nr |
    head -n 10
```

### Download sequences and regenerate lineage information

```bash
mkdir -p ~/data/bacteria/GENOMES
cd ~/data/bacteria/GENOMES

cat ../summary/ABBR.csv |
    grep -v '^#' |
    perl -nla -F"," -e 'print qq{$F[0],$F[7]}' |
    uniq |
    perl ~/Scripts/withncbi/taxon/strain_info.pl --stdin --withname --file bac_ncbi.csv

echo "#strain_name,accession,strain_taxon_id" > bac_name_acc_id.csv
cat ../summary/ABBR.csv |
	grep -v '^#' |
    perl -nla -F"," -e '
        my $acc = $F[6];
        $acc =~ s/"//g;
        $acc =~ s/\.\d+//g;
        for my $s (split /\|/, $acc) {
            print qq{$F[7],$s,$F[0]};
        }
    ' |
    sort >> bac_name_acc_id.csv
cat bac_name_acc_id.csv | wc -l
```

```bash
cat bac_name_acc_id.csv |
    grep -v '^#' |
    2>&1 parallel --colsep ',' --no-run-if-empty --linebuffer -k -j 3 "
        echo -e '==> id: [{1}]\tseq: [{2}]\n'
        mkdir -p {1}
        if [[ -e '{1}/{2}.gff' && -e '{1}/{2}.fa' ]] ; then
            echo -e '    Sequence [{1}/{2}] exists, next\n'
            exit
        fi

        # gb
        echo -e '    [{1}/{2}].gb'
        curl -Ls \
            'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={2}&rettype=gb&retmode=text' \
            > {1}/{2}.gb

        # fasta
        echo -e '    [{1}/{2}].fa'
        curl -Ls \
            'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={2}&rettype=fasta&retmode=text' \
            > {1}/{2}.fa

        # gff
        echo -e '    [{1}/{2}].gff'
        perl ~/Scripts/withncbi/taxon/bp_genbank2gff3.pl {1}/{2}.gb -o stdout > {1}/{2}.gff
        perl -i -nlp -e '/^\#\#FASTA/ and last' {1}/{2}.gff

        echo
    " |
    tee bac_seq.log
    
# count downloaded sequences
find . -maxdepth 2 -name "*.fa" | wc -l
find . -maxdepth 1 -type d | wc -l

# failed files
find . -maxdepth 2 -type f -size -1k | grep ".fa$"
```

### Numbers for different ranks -- Genus and species

```bash
cd ~/data/bacteria/summary/

# count every ranks
cat ABBR.csv | sed -e '1d' | cut -d',' -f 3 | sort | uniq > species.list
cat ABBR.csv | sed -e '1d' | cut -d',' -f 4 | sort | uniq > genus.list.tmp
cat ABBR.csv | sed -e '1d' | cut -d',' -f 5 | sort | uniq > subgroup.list
wc -l subgroup.list genus.list.tmp species.list
# because we only needs Staphylococcus, so in subgroup and genus, we could only find 1, and in species we got 11.

rm *.tmp
```



## Genomes from bacteria_ar: assembly

```bash
export RANK_NAME=bacteria_ar

mkdir -p ~/data/alignment/${RANK_NAME}        # Working directory
cd ~/data/alignment/${RANK_NAME}

mysql -ualignDB -palignDB ar_genbank -e "
    SELECT
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND taxonomy_id != species_id               # no strain ID
        AND genus_id = 1279
    " \
    > raw.tsv

cat raw.tsv |
    grep -v '^#' |
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "1,2,3" -s '\t' -m 3 --shortsub |
    (echo -e '#name\tftp_path\torganism\tassembly_level' && cat ) |
    perl -nl -a -F"," -e '
        BEGIN{my %seen};
        /^#/ and print and next;
        /^organism_name/i and next;
        $seen{$F[5]}++;
        $seen{$F[5]} > 1 and next;
        printf qq{%s\t%s\t%s\t%s\n}, $F[5], $F[3], $F[1], $F[4];
        ' |
    keep-header -- sort -k3,3 -k1,1 \
    > ${RANK_NAME}.assembly.tsv

# comment out unneeded assembly levels

# find potential duplicated strains or assemblies
cat ${RANK_NAME}.assembly.tsv |
    cut -f 1 |
    sort |
    uniq -c |
    sort -nr

# Cleaning
rm raw*.*sv

unset RANK_NAME
```

```bash
perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
	-f bacteria_ar.assembly.tsv \
	-o ASSEMBLY

bash ASSEMBLY/bacteria_ar.assembly.rsync.sh
# check whether bacteria_ar downloaded successfully
bash ASSEMBLY/bacteria_ar.assembly.collect.sh
```

```bash
cd ~/data/alignment/bacteria_ar
# count numbers of different species in genus Staphylococcus
# grep -Ev could filter the blank and comment lines
cat ASSEMBLY/bacteria_ar.assembly.collect.csv | cut -d "," -f 2 | grep -v "Organism_name" | grep -Ev "^$|[#;]" | cut -d " " -f 1,2 | sort | uniq -c
```



```bash
cat bacteria_ar.assembly.tsv | cut -f 3 | sed '1d' | uniq -c > bacteria_ar_info.tsv
# count different Staphylococcus genomes
```

### move genomic to DB

```bash
cd ~/data/alignment/bacteria_ar/ASSEMBLY
for FILE_NAME in `ls`
do
    cd ${FILE_NAME}
    cp *genomic.fna.gz ~/data/alignment/bacteria_ar/DB2
    cd ~/data/alignment/bacteria_ar/ASSEMBLY
done
echo "Copy completed"
```



#  Build db of genomes for blast

##  Generate database of Staphylococcus in blast

```bash
cd ~/data/bacteria/gr_db

makeblastdb -in Sta_genus.fa -dbtype nucl -parse_seqids \
-out Sta_genus -title "Staphylococcus_Genome_DataBase"
```



##  Align primer within database using blastn

```bash
cd ~/data/bacteria/gr_db

blastn -query ../biofire/S_aureus_primer.fasta -db Sta_genus -task blastn -out test_btop.txt -num_threads 6
blastn -query ../biofire/S_epidermidis_primer.fasta -db Sta_genus -task blastn -outfmt 7 -out ../biofire/test.txt -num_threads 6
```



##  Decompress all the genomic fasta to one file

```bash
gunzip -c -d *.gz >> ../Sta_genus_db/db.fna
```

```bash
makeblastdb -in db.fna -dbtype nucl -parse_seqids \
-out Sta_genus_ar -title "StaGenomeDatabase"
```

```bash
cd ~/data/alignment/bacteria_ar/Sta_genus_db

blastn -query ~/data/biofire/Sau-rpoB-F.fa -db Sta_genus_ar \
-task blastn -outfmt 7 -max_target_seqs 10000 \
-out ~/data/biofire/blastn/Sau-rpoB-F.tsv

blastn -query ~/data/biofire/Sau-rpoB-R.fa -db Sta_genus_ar \
-task blastn -outfmt 7 -max_target_seqs 10000 \
-out ~/data/biofire/blastn/Sau-rpoB-R.tsv

blastn -query ~/data/biofire/Sau-gyrB-F.fa -db Sta_genus_ar \
-task blastn -outfmt 7 -max_target_seqs 10000 \
-out ~/data/biofire/blastn/Sau-gyrB-F.tsv

blastn -query ~/data/biofire/Sau-gyrB-R.fa -db Sta_genus_ar \
-task blastn -outfmt 7 -max_target_seqs 10000 \
-out ~/data/biofire/blastn/Sau-gyrB-R.tsv

blastn -query ~/data/biofire/Sau-amplify-F.fa -db Sta_genus_ar \
-task blastn -outfmt 7 -max_target_seqs 10000 \
-out ~/data/biofire/blastn/Sau-amplify-F.tsv

blastn -query ~/data/biofire/Sau-amplify-R.fa -db Sta_genus_ar \
-task blastn -outfmt 7 -num_descriptions 10000 \
-out ~/data/biofire/blastn/Sau-amplify-R.tsv
```



```bash
cat Sau-gyrB-R.tsv | sed "1,5d" | sed "1i\#qacc\tsacc\tidentity\talilength\tmismatches\tgap\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > Sau-gyrB-R1.tmp.tsv

tsv-filter -H --eq identity:100.000 --eq alilength:24 --eq mismatches:0 Sau-gyrB-R1.tsv > Sau-gyrB-R2.tmp.tsv
```

```bash
sed -i '1i\sacc\tstrains' bac_acc_info.tsv

tsv-select Sau-gyrB-F2.tmp.tsv Sau-gyrB-R2.tmp.tsv -H -f 2 > Sau-gyrB.tmp.tsv
tsv-uniq -H --a 2 Sau-gyrB.tmp.tsv > Sau-gyrB-1.tmp.tsv
tsv-join -H -f bac_acc_info.tsv -k sacc --append-fields strains Sau-gyrB-1.tmp.tsv > Sau-gyrB-info.tsv
cat Sau-gyrB-info.tsv | cut -f 2 | sed '1d' | sort | uniq -c > Sau-gyrB.tsv
```



```bash
mysql -ualignDB -palignDB ar_refseq -e "
    SELECT
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND species_id in (985002, 1654388)
    " \
    > raw.tsv
    
mysql -ualignDB -palignDB ar_genbank -e "
    SELECT
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND taxonomy_id != species_id               # no strain ID
        AND species_id in (985002, 1654388)
    " \
    >> raw.tsv
```

```bash
perl ../bin/blast+/update_blastdb.pl --passive --decompress ref_prok_rep_genomes
# --passive: Use passive FTP, useful when behind a firewall or working in the cloud (default:true)
```

```bash
cd /mnt/e/db
blastn -query ~/data/biofire/Pae-amplify-F.fa -db ref_prok_rep_genomes \
-task blastn -outfmt 7 -max_target_seqs 10000 \
-out /mnt/c/Users/59717/Desktop/Blast/Pae-amp-F.tsv
```





```bash
samtools faidx Pseudomonas_aeruginosa.fasta NC_002516.2:1116019-1116038 > Pae-amp.fasta
blastn -query ../Pae-amp.fasta -db Genomes -task blastn -out ../Pae-amp.txt -num_threads 6
```

