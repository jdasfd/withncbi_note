# Aligning HIV species-specific DNA-tags

## NCBI Assembly

```bash
mkdir -p /mnt/e/data/alignment/Viruses
cd /mnt/e/data/alignment/Viruses

mysql -ualignDB -palignDB ar_refseq -e "
    SELECT 
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar 
    WHERE 1=1
        AND taxonomy_id != species_id               # no strain ID
        AND genus_id = 11646
    " \
    > raw.tsv


mysql -ualignDB -palignDB ar_genbank -e "
    SELECT 
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar 
    WHERE 1=1
        AND taxonomy_id != species_id               # no strain ID
        AND genus_id = 11646
    " \
    >> raw.tsv

cat raw.tsv |
    grep -v '^#' |
    perl /mnt/e/Scripts/withncbi/taxon/abbr_name.pl -c "1,2,3" -s '\t' -m 3 --shortsub |
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
    > Viruses.assembly.tsv

cat Viruses.assembly.tsv | cut -f 1 | sed -e '1d' | sort | uniq -c | sort -nr

rm raw*.*sv
```

```bash
cd /mnt/e/data/alignment/Viruses

cp Viruses.assembly.tsv /mnt/e/Scripts/withncbi/pop

perl /mnt/e/Scripts/withncbi/taxon/assembly_prep.pl \
    -f /mnt/e/Scripts/withncbi/pop/Viruses.assembly.tsv \
    -o ASSEMBLY

bash ASSEMBLY/Viruses.assembly.rsync.sh
bash ASSEMBLY/Viruses.assembly.collect.sh
```

## Count strains

```bash
cd /mnt/e/data/alignment/Viruses

for dir in $(find ASSEMBLY -maxdepth 1 -mindepth 1 -type d | sort); do
    1>&2 echo "==> ${dir}"
    name=$(basename ${dir})
    
    find ${dir} -type f -name "*_genomic.fna.gz" |
        grep -v "_from_" |
        xargs cat |
        faops n50 -C -S stdin |
        (echo -e "name\t${name}" && cat) |
        datamash transpose
done |
    tsv-uniq |
    tee ASSEMBLY/n50.tsv

cat ASSEMBLY/n50.tsv |
    tsv-filter \
        -H --or \
        --le 4:2000 \
        --ge 2:100000 |
    tr "\t" "," \
    > ASSEMBLY/n50.pass.csv

wc -l ASSEMBLY/n50*

tsv-join \
    ASSEMBLY/Viruses.assembly.collect.csv \
    --delimiter "," -H --key-fields 1 \
    --filter-file ASSEMBLY/n50.pass.csv \
    > ASSEMBLY/Viruses.assembly.pass.csv

wc -l ASSEMBLY/Viruses.assembly*csv
```

```bash
cd /mnt/e/data/alignment/Viruses

parallel --no-run-if-empty --linebuffer -k -j 4 '
    n_species=$(cat ASSEMBLY/Viruses.assembly.pass.csv |
        cut -d"," -f 2 |
        grep -v "Candidatus" |
        grep "{}" |
        cut -d" " -f 1,2 |
        sort |
        uniq |
        wc -l)
    
    n_strains=$(cat ASSEMBLY/Viruses.assembly.pass.csv |
        cut -d"," -f 2 |
        grep -v "Candidatus" |
        grep "{}" |
        cut -d" " -f 1,2 |
        sort |
        wc -l)
    
    printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    ' ::: $(
        cat ASSEMBLY/Viruses.assembly.pass.csv |
            sed -e '1d' |
            cut -d"," -f 2 |
            grep -v "Candidatus" |
            cut -d" " -f 1 |
            sort |
            uniq
    )
```

## Raw phylogenetic tree by MinHash

```bash
mkdir -p /mnt/e/data/alignment/Viruses/mash
cd /mnt/e/data/alignment/Viruses/mash

for name in $(cat ../ASSEMBLY/Viruses.assembly.pass.csv | sed -e '1d' | cut -d"," -f 1 ); do
    2>&1 echo "==> ${name}"
    
    if [[ -e ${name}.msh ]]; then
        continue
    fi
    
    find ../ASSEMBLY/${name} -name "*.fsa_nt.gz" -or -name "*_genomic.fna.gz" |
        grep -v "_from_" |
        xargs cat |
        mash sketch -k 21 -s 100000 -p 10 - -I "${name}" -o ${name}
done

mash triangle -E -p 10 -l <(
    cat ../ASSEMBLY/Viruses.assembly.pass.csv | sed -e '1d' | cut -d"," -f 1 | parallel echo "{}.msh"
    ) \
    > dist.tsv

tsv-select -f 1-3 dist.tsv |
    (tsv-select -f 2,1,3 dist.tsv && cat) |
    (
        cut -f 1 dist.tsv |
            tsv-uniq |
            parallel -j 1 --keep-order 'echo -e "{}\t{}\t0"' &&
        cat
    ) \
    > dist_full.tsv

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
        
        group <- cutree(clusters, h=0.5) # k=3
        groups <- as.data.frame(group)
        groups$ids <- rownames(groups)
        rownames(groups) <- NULL
        groups <- groups[order(groups$group), ]
        write_tsv(groups, "groups.tsv")
    '

nw_display -s -b 'visibility:hidden' -w 600 -v 30 tree.nwk |
    rsvg-convert -o /mnt/e/data/alignment/Viruses/mash/Viruses.png
```

## Groups and targets

We chose a  genome first sequenced as a reference.

```bash
mkdir /mnt/e/data/alignment/Viruses/taxon
cd /mnt/e/data/alignment/Viruses/taxon

cp ../mash/tree.nwk ../mash/groups.tsv .

ARRAY=('HIV::L_Hum_HIV_1_CRF03_AB')

echo -e "#Serial\tGroup\tCount\tTarget\tSequencing" > group_target.tsv

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    TARGET_NAME="${item##*::}"
    
    SERIAL=$(
        cat ../mash/groups.tsv |
            tsv-filter --str-eq 2:${TARGET_NAME} |
            tsv-select -f 1
    )

    cat ../mash/groups.tsv |
        tsv-filter --str-eq 1:${SERIAL} |
        tsv-select -f 2 \
        > ${GROUP_NAME}

    COUNT=$(cat ${GROUP_NAME} | wc -l )

    echo -e "${SERIAL}\t${GROUP_NAME}\t${COUNT}\t${TARGET_NAME}\t" >> group_target.tsv

done

mlr --itsv --omd cat group_target.tsv

cat <<'EOF' > chr-level.list
L_Hum_HIV_1_CRF03_AB
EOF
```

## HIV: prepare

```bash
cd /mnt/e/data/alignment/Viruses

# prep
egaz template \
    ASSEMBLY \
    --prep -o GENOMES \
    $( cat taxon/group_target.tsv | sed -e '1d' | cut -f 4 | parallel -j 1 echo " --perseq {} " ) \
    $( cat taxon/chr-level.list | parallel -j 1 echo " --perseq {} " ) \
    --min 5000 --about 5000000 \
    -v --repeatmasker "--species Viruses --parallel 12"

bash GENOMES/0_prep.sh

# gff
for n in $(cat taxon/group_target.tsv | sed -e '1d' | cut -f 4 ) \
    $( cat taxon/chr-level.list ) \
    ; do
    FILE_GFF=$(find ASSEMBLY -type f -name "*_genomic.gff.gz" | grep "${n}")
    echo >&2 "==> Processing ${n}/${FILE_GFF}"
    
    gzip -d -c ${FILE_GFF} > GENOMES/${n}/chr.gff
done
```

## HIV: run

```bash
cd /mnt/e/data/alignment/Viruses/

# sanger
egaz template \
    GENOMES/L_Hum_HIV_1_CRF03_AB \
    GENOMES/L_Hum_HIV_1_CRF04_cpx \
    --multi -o groups/HIV/ \
    --multiname sanger \
    --tree taxon/tree.nwk \
    --parallel 12 -v

bash groups/HIV/1_pair.sh
bash groups/HIV/3_multi.sh

# multi
egaz template \
    GENOMES/L_Hum_HIV_1_CRF03_AB \
    $(cat taxon/HIV | grep -v -x "L_Hum_HIV_1_CRF03_AB" | parallel -j 1 echo "GENOMES/{}") \
    --multi -o groups/HIV/ \
    --tree taxon/tree.nwk \
    --parallel 12 -v

bash groups/HIV/1_pair.sh
bash groups/HIV/3_multi.sh
```

## fas from axt files

```bash
cd /mnt/e/data/alignment/Viruses/groups/HIV/Pairwise

foreach dir in $(find -maxdepth 1 -mindepth 1 -type d); do
    name=$(basename ${dir})
    p=${name%%vs*}
    q=${name##vs*}
    fasops axt2fas ./${name}/axtNet/*.axt.gz -o ./${name}/axtNet/result.fasta \
        -l 1000 -t ${p} -q ${q} \
        -s /mnt/e/data/alignment/Viruses/GENOMES/${q}/chr.sizes
done
```
