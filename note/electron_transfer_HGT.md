# Electron transfer related proteins

```bash
cat ASSEMBLY/Pseudomonas.assembly.pass.csv |
  tsv-summarize -H -d ',' -g RefSeq_category --count
#RefSeq_category,count
#Representative Genome,533
#,1404
#Reference Genome,15
```

## Background info

### Software

- Install `nwr` and create a local taxonomy and assembly database.

```bash
brew install wang-q/tap/nwr # 0.5.5 or above
brew install sqlite         # 3.34 or above

nwr download
nwr txdb

nwr ardb
nwr ardb --genbank
```

`hmmsearch`: Search a protein profile HMM against a protein sequence database.

`hmmscan`: Search a protein sequence against a protein profile HMM database.

That is the difference of the two methods.

## Proteins info

The two proteins were selected from my teacher's protocol in [Pseudomonas.md](Pseudomonas.md). The link is directed to my study records. You can find the table [here](Pseudomonas.md#interproscan-on-all-proteins-of-typical-strains)

- IPR000813: 7Fe ferredoxin

  7Fe ferredoxin are iron-sulphur proteins (铁硫蛋白) that mediate electron transfer in a range of metabolic reactions. GO:0009055

- IPR024922: Rubredoxin

  Rubredoxin is a low molecular weight iron-containing bacterial and archaeal protein that is involved in electron transfer, sometimes replacing ferredoxin as an electron carrier. GO:0009055, GO:0005506

## IPR024922 - Rubredoxin

- PIRSF: PIRSF000071 Rubredoxin
- Pfam: PF00301: Rubredoxin
- Panther: PTHR47627: Rubredoxin

In PIRSF database, the same Pfam domain was contained, so I used the Pfam HMM.

### hmmsearch for Rubredoxin

- `hmmsearch` for Rubredoxin domains in each strain

The goal of the step is to search target protein from each strain according to HMM.

```bash
cd /mnt/e/data/Pseudomonas

mkdir -p Rubr/HMM

curl -L https://pfam.xfam.org/family/PF00301/hmm > Rubr/HMM/Rubr.hmm
#curl -L https://www.pantherdb.org/panther/exportHmm.jsp?acc=PTHR47627 > Rubr/HMM/PTHR47627.hmm

# using Rubredoxin.hmm to search in each strains
E_VALUE=1e-20
for domain in Rubr; do
    >&2 echo "==> domain [${domain}]"

    if [ -e Rubr/${domain}.replace.tsv ]; then
        continue;
    fi

    for ORDER in $(cat order.lst); do
        >&2 echo "==> ORDER [${ORDER}]"

        cat taxon/${ORDER} |
            parallel --no-run-if-empty --linebuffer -k -j 8 "
                gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                    hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw Rubr/HMM/${domain}.hmm - |
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
        > Rubr/${domain}.replace.tsv

    >&2 echo
done

# 2 protein replace file intersection
#tsv-join Rubr/Rubr.replace.tsv \
#    -f Rubr/PTHR47627.replace.tsv \
#    > Rubr/Rubredoxin.replace.tsv

wc -l Rubr/Rubredoxin.replace.tsv
#2089 Rubr/Rubredoxin.replace.tsv
# remember that 1952 strains passed were obtained

# extract all Rubr pro_seqs
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 Rubr/Rubredoxin.replace.tsv) Rubr/Rubr.fa

muscle -in Rubr/Rubr.fa -out Rubr/Rubr.aln.fa

FastTree Rubr/Rubr.aln.fa > Rubr/Rubr.aln.newick

nw_reroot Rubr/Rubr.aln.newick $(nw_labels Rubr/Rubr.aln.newick | grep -E "Bac_subti|Sta_aure") |
    nw_order -c n - \
    > Rubr/Rubr.reroot.newick

# Rubredoxin among all reference genomes
cat Rubr/Rubredoxin.replace.tsv | grep -v 'GCF'
#YP_004994544.1  Acin_pittii_PHEA_2_YP_004994544
#NP_820858.1     Co_burn_RSA_493_NP_820858
#NP_417190.1     Es_coli_K_12_MG1655_NP_417190
#NP_311593.1     Es_coli_O157_H7_Sakai_NP_311593
#YP_005228417.1  Kle_pneumon_pneumoniae_HS11286_YP_005228417
#NP_461761.1     Salm_enterica_enterica_Typhimurium_LT2_NP_461761
#NP_708517.1     Shig_fle_2a_301_NP_708517
#NP_254038.1     Pseudom_aeru_PAO1_NP_254038
#NP_254037.1     Pseudom_aeru_PAO1_NP_254037
#YP_002517953.1  Cau_vib_NA1000_YP_002517953
#NP_217767.1     My_tube_H37Rv_NP_217767
#NP_217768.1     My_tube_H37Rv_NP_217768
```

- Counting copies of each strain

```bash
cd /mnt/e/data/Pseudomonas

cat PROTEINS/all.strain.tsv |
    sed '1d' |
    grep -F -f <(cut -f 2 Rubr/Rubredoxin.replace.tsv) |
    cut -f 2 |
    tsv-summarize -g 1 --count \
    > Rubr/strains.copy.tsv

wc -l Rubr/strains.copy.tsv
#1567 Rubr/strains.copy.tsv
# all strains passed up to 1952, so there are strains missing the Rubredoxin

cat strains.lst |
    grep -v -F -f <(cat Rubr/strains.copy.tsv | cut -f 1) |
    awk '{print $0"\t"0}' \
    >> Rubr/strains.copy.tsv

wc -l Rubr/strains.copy.tsv
#1952 Rubr/strains.copy.tsv

(echo -e "strains\tcopy_num" && cat Rubr/strains.copy.tsv) > \
    temp && mv temp Rubr/strains.copy.tsv

cat strains.taxon.tsv |
    cut -f 1,4 |
    sed '1istrains\tspecies' |
    tsv-join -H --filter-file Rubr/strains.copy.tsv -k strains --append-fields copy_num \
    > Rubr/species.copy.tsv

# species with Rubredoxin more than 1 copy
cat Rubr/species.copy.tsv |
    tsv-filter -H --ge copy_num:2 |
    tsv-select -f 2 |
    tsv-summarize -H -g species --count \
    > Rubr/species.2copy.tsv

cat ASSEMBLY/Pseudomonas.assembly.pass.csv |
    sed -e '1d' |
    tsv-select -d, -f 3 |
    nwr append stdin -r species |
    tsv-summarize -g 2 --count \
    > Rubr/species.count.tsv

# all strains are 2 copy in a species
tsv-join --filter-file Rubr/species.count.tsv \
-k 1 --append-fields 2 Rubr/species.2copy.tsv |
    tsv-filter --ff-eq 2:3 |
    mlr --itsv --omd cat
```

| Acidihalobacter aeolianus         | 1   | 1_2 |
|-----------------------------------|-----|-----|
| Acidihalobacter yilgarnensis      | 1   | 1   |
| Alkalilimnicola ehrlichii         | 1   | 1   |
| Azotobacter chroococcum           | 3   | 3   |
| Legionella longbeachae            | 1   | 1   |
| Legionella sainthelensi           | 1   | 1   |
| Marinobacter nauticus             | 2   | 2   |
| Methylocaldum marinum             | 1   | 1   |
| Methylogaea oryzae                | 1   | 1   |
| Methylomicrobium album            | 1   | 1   |
| Methylomonas denitrificans        | 1   | 1   |
| Methylomonas koyamae              | 1   | 1   |
| Methylomonas methanica            | 1   | 1   |
| Methylomonas rhizoryzae           | 1   | 1   |
| Methylotuvimicrobium alcaliphilum | 1   | 1   |
| Methylotuvimicrobium buryatense   | 1   | 1   |
| Methylovulum psychrotolerans      | 1   | 1   |
| Moraxella catarrhalis             | 16  | 16  |
| Mycobacterium tuberculosis        | 1   | 1   |
| Plesiomonas shigelloides          | 1   | 1   |
| Pseudomonas aeruginosa            | 391 | 391 |
| Pseudomonas alcaligenes           | 3   | 3   |
| Pseudomonas citronellolis         | 2   | 2   |
| Pseudomonas knackmussii           | 2   | 2   |
| Pseudomonas lalkuanensis          | 2   | 2   |
| Pseudomonas otitidis              | 3   | 3   |
| Pseudomonas sessilinigenes        | 2   | 2   |
| Pseudoxanthomonas spadix          | 1   | 1   |
| Spiribacter curvatus              | 1   | 1   |
| Sulfurivermis fontis              | 1   | 1   |
| Thioalkalivibrio sulfidiphilus    | 1   | 1   |

### Compare protein seqs with PFAM using hmmscan

TIGRFAM/PFAM database could be seen in [hmm.md](hmm.md)

The previous step: [hmmsearch for Rubredoxin](#hmmsearch-for-rubredoxin) used the HMMs of Rubredoxin domain to search against strain protein files and finally got all proteins with that domain. This step is adopted to exclude those proteins with the domain but not our targets.

The reason of using PFAM database: only take into account of domain HMM from PFAM db.

- `hmmscan` for strains filtering

```bash
cd /mnt/e/data/Pseudomonas

cat ~/data/HMM/TIGRFAM/HMM/*.HMM > Rubr/HMM/TIGRFAM.hmm
cat ~/data/HMM/PGAP/HMM/hmm_PGAP/*.HMM > Rubr/HMM/PGAP.hmm
cp ~/data/HMM/PFAM/*.hmm ~/data/HMM/PFAM/*.hmm.dat ~/data/HMM/PFAM/*.hmm.h* Rubr/HMM/

# hmmpress to prepare an database
hmmpress Rubr/HMM/TIGRFAM.hmm
hmmpress Rubr/HMM/PGAP.hmm
# hmmpress could prepare an HMM database for faster hmmscan searches

# faops extract protein seq
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 Rubr/Rubredoxin.replace.tsv) Rubr/Rubr.fa

# using protein seqs to scan directly in PFAM database
E_VALUE=1e-10
#if [ ! -s Rubr/HMM/PGAP.hmm ]; then
#    echo no PGAP
#if [ ! -s Rubr/HMM/PGAP.hmm ]; then
#    echo no TIGRFAM
if [ ! -s Rubr/HMM/Pfam-A.hmm ]; then
    echo no Pfam
    exit
else
    echo >&2 hmmscan start: Rubr
    hmmscan --cpu 6 -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw \
        -o Rubr/Rubr.pfam.txt --tblout Rubr/Rubr.pfam.tbl \
        Rubr/HMM/Pfam-A.hmm Rubr/Rubr.fa
    echo hmmscan complete
fi
#    hmmscan --cpu 6 -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw \
#        -o Rubr/Rubr.pgap.txt --tblout Rubr/Rubr.pgap.tbl \
#        Rubr/HMM/PGAP.hmm Rubr/Rubr.fa
#    echo hmmscan complete
#    hmmscan --cpu 6 -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw \
#        -o Rubr/Rubr.tigrfam.txt --tblout Rubr/Rubr.tigrfam.tbl \
#        Rubr/HMM/TIGRFAM.hmm Rubr/Rubr.fa
#    echo hmmscan complete

# reformat and extract results from tbl
cat Rubr/Rubr.pfam.tbl |
    grep -v '^#' |
    perl -nl -e '$string = $_;@p = ();
    for ($i = 1; $i <= 19; $i++){
    if($i != 19){
        $string=~s/^(.+?)\s+//;
        push(@p, $1);
    }
    else{
        push(@p, $string);
    }}
    print join ("\t", @p);
    ' > Rubr/Rubr.pfam.tsv

cat Rubr/Rubr.pfam.tsv |
    tsv-select -f 1 |
    tsv-uniq
#Rubredoxin
#Pyr_redox_2
#Pyr_redox
#Pyr_redox_3
#Rbx_binding
#Lycopene_cycl
# This means hmmscan will give you more than 1 domains

cat Rubr/Rubr.fa | grep '^>' | wc -l
#1990
# Total 1990 protein seqs were used for hmmscan

cat Rubr/Rubr.pfam.tsv | grep -v '^#' | wc -l
#2027
# Got 2027 lines, one protein may had more than 1 hit by hmmscan

cat Rubr/Rubr.pfam.tsv |
    tsv-select -f 3 |
    tsv-uniq |
    wc -l
#1990
# All proteins completed hmmscan

cat Rubr/Rubr.pfam.tsv | tsv-select -f 3 | tsv-uniq --repeated
#Alk_ehr_MLHE_1_GCF_000014785_1_WP_011628506
#Alt_mac_GCF_903772925_1_WP_179983007
#Halot_nea_c2_GCF_000024765_1_WP_012825058
#Mora_cat_GCF_000766665_1_WP_003666012
#Mora_cat_GCF_002073215_2_WP_063453878
#Mora_cat_GCF_003443915_1_WP_003667335
#Mora_cat_GCF_003443975_1_WP_117905972
#Mora_cat_GCF_003971285_1_WP_003666012
#Mora_cat_GCF_003971345_1_WP_003666012
#Mora_cat_GCF_003971365_1_WP_003666012
#Mora_cat_GCF_900476075_1_WP_063453878
#Mora_osl_GCF_002753715_1_WP_100271360
#Thioa_versu_GCF_001020955_1_WP_047251304
# with more than 1 hit after hmmscan

cat Rubr/Rubr.pfam.tsv |
    tsv-summarize -g 1 --count |
    sort -r -nk 2,2
#Rubredoxin      1990
#Pyr_redox_2     13
#Pyr_redox       13
#Pyr_redox_3     9
#Rbx_binding     1
#Lycopene_cycl   1
# Got 1990 Rubredoxin = Rubr.fa num of seqs

# keep the minimum E-value strain-protein id and move them into a LIST
cat Rubr/Rubr.pfam.tsv |
    tsv-select -f 3,5 |
    tsv-summarize -g 1 --min 2 |
    cut -f 1 \
    > Rubr/Rubr.Emin.lst

wc -l Rubr/Rubr.Emin.lst
#1990 Rubr/Rubr.Emin.lst
```

- Counting copies of each strain

```bash
cat PROTEINS/all.strain.tsv |
    sed '1d' |
    grep -F -f <(cut -f 2 Rubr/Rubredoxin.replace.tsv) |
    tsv-join --filter-file Rubr/Rubr.Emin.lst -k 1 |
    cut -f 2 |
    tsv-summarize -g 1 --count \
    > Rubr/strains.Emin.copy.tsv

wc -l Rubr/strains.Emin.copy.tsv
#1502 Rubr/strains.Emin.copy.tsv
# all strains passed up to 1952, so there are strains missing the Rubredoxin

cat strains.lst |
    grep -v -F -f <(cat Rubr/strains.Emin.copy.tsv | cut -f 1) |
    awk '{print $0"\t"0}' \
    >> Rubr/strains.Emin.copy.tsv

wc -l Rubr/strains.Emin.copy.tsv
#1952 Rubr/strains.Emin.copy.tsv

(echo -e "strains\tcopy_num" && cat Rubr/strains.Emin.copy.tsv) > \
    temp && mv temp Rubr/strains.Emin.copy.tsv

cat strains.taxon.tsv |
    cut -f 1,4 |
    sed '1istrains\tspecies' |
    tsv-join -H --filter-file Rubr/strains.Emin.copy.tsv -k strains --append-fields copy_num \
    > Rubr/species.Emin.copy.tsv

# species with Rubredoxin more than 1 copy
cat Rubr/species.Emin.copy.tsv |
    tsv-filter -H --ge copy_num:2 |
    tsv-select -f 2 |
    tsv-summarize -H -g species --count \
    > Rubr/species.Emin.2copy.tsv

cat ASSEMBLY/Pseudomonas.assembly.pass.csv |
    sed -e '1d' |
    tsv-select -d, -f 3 |
    nwr append stdin -r species |
    tsv-summarize -g 2 --count \
    > Rubr/species.Emin.count.tsv

# all strains are 2 copy in a species
tsv-join --filter-file Rubr/species.Emin.count.tsv \
-k 1 --append-fields 2 Rubr/species.Emin.2copy.tsv |
    tsv-filter --ff-eq 2:3
#Acidihalobacter yilgarnensis    1       1
#Alkalilimnicola ehrlichii       1       1
#Legionella longbeachae  1       1
#Legionella sainthelensi 1       1
#Methylocaldum marinum   1       1
#Methylogaea oryzae      1       1
#ethylomicrobium album  1       1
#Methylomonas denitrificans      1       1
#Methylomonas koyamae    1       1
#Methylomonas methanica  1       1
#Methylotuvimicrobium alcaliphilum       1       1
#Methylotuvimicrobium buryatense 1       1
#Pseudomonas aeruginosa  391     391
#Pseudomonas alcaligenes 3       3
#Pseudomonas citronellolis       2       2
#Pseudomonas knackmussii 2       2
#Pseudomonas lalkuanensis        2       2
#Pseudomonas otitidis    3       3
#Pseudomonas sessilinigenes      2       2
#Pseudoxanthomonas spadix        1       1
#Spiribacter curvatus    1       1
#Sulfurivermis fontis    1       1
#Thioalkalivibrio sulfidiphilus  1       1
```

### Multiple rounds using `blastp`

The step is adopted for blasting all proteins (filtered by `hmmscan`) against all protein seqs again to extract those proteins missed during `hmmsearch`. The goal is to acquire proteins directly from sequences (not just relied on `hmmsearch`).

```bash
cd /mnt/e/data/Pseudomonas
mkdir -p Rubr/blastp

# remove null character
cat PROTEINS/all.replace.fa |
    sed 's/\x0//g' \
    > Rubr/blastp/all.db.fa

# extract seqs with minimum E-value
faops some PROTEINS/all.replace.fa Rubr/Rubr.Emin.lst Rubr/blastp/Rubr.Emin.fa
faops size Rubr/blastp/Rubr.Emin.fa | wc -l
#1983

# makeblastdb - all pro_seqs
makeblastdb -in Rubr/blastp/all.db.fa -dbtype prot -out Rubr/blastp/all_pro

# blastp 1st
E_VALUE=1e-20
blastp -db Rubr/blastp/all_pro -query Rubr/blastp/Rubr.Emin.fa \
-out Rubr/blastp/Rubr.1st.tsv -outfmt 6 -evalue ${E_VALUE} -num_threads 6

cat Rubr/blastp/Rubr.1st.tsv |
    cut -f 2 |
    sort -n |
    uniq |
    wc -l
#3260
# after blastp, we scanned more proteins

# extract seqs from 1st blastp result
faops some PROTEINS/all.replace.fa <(cut -f 2 Rubr/blastp/Rubr.1st.tsv | sort -n | uniq) Rubr/blastp/Rubr.1st.fa
faops size Rubr/blastp/Rubr.1st.fa | wc -l
#
```

- Build tree by `iTOL` online

The previous step provide us the list of all more than 1 copy strains among the species. According to the purpose, the protein tree and species tree are built and are compared to find out whether the topological structure of two trees are different.

```bash
cd /mnt/e/data/Pseudomonas

# provide iTOL an TXT annotation file to change color
# genus of Pseudomonas will be colored by red
cat strains.taxon.tsv |
    tsv-select -f 1,4 |
    tsv-filter --str-in-fld 2:Pseudomonas |
    tsv-select -f 1 |
    awk '{print $0"\tlabel\t#f44336"}' \
    > Rubr/tree/strains.label.txt

(echo -e "TREE_COLORS\nSPEARATOR TAB\nDATA" && cat Rubr/tree/strains.label.txt) \
    > temp && mv temp Rubr/tree/strains.label.txt
```

## IPR000813 - 7Fe ferredoxin

- PRINTS: PR00354 7Fe ferredoxin signature

PRINTS database does not have the HMM of the protein.

7FE8SFRDOXIN is a 3-element fingerprint that provides a signature for 7Fe ferredoxins. The fingerprint was derived from an initial alignment of 2 sequences: motif 1 contains 2 of 3 invariant Cys residues contributing to the 3Fe-4S cluster, which is located upstream from the 4Fe-4S domain; motif 2 contains the first invariant Cys of the 4Fe-4S cluster; and motif 3 contains 3 invariant cysteines of the 4Fe-4S cluster and one Cys ligand of its downstream 3Fe-4S counterpart (see scheme above).

- Pfam: PF00037: Fer4, which is the 4Fe-4S motif 1

```bash
cat STRAINS/Pseudom_aeru_PAO1/*.tsv |
    grep "IPR000813"
# this step will provide you all the 7Fe ferredoxin info in the previous analyze
# the 7Fe ferredoxin was recorded in the PRINTS database

mkdir -p Ferr/HMM

curl -L https://pfam.xfam.org/family/PF00037/hmm > Ferr/HMM/Fer4.hmm

E_VALUE=1e-20
for domain in Fer4; do
    >&2 echo "==> domain [${domain}]"

    if [ -e Rubr/${domain}.replace.tsv ]; then
        continue;
    fi

    for ORDER in $(cat order.lst); do
        >&2 echo "==> ORDER [${ORDER}]"

        cat taxon/${ORDER} |
            parallel --no-run-if-empty --linebuffer -k -j 8 "
                gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                    hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw Ferr/HMM/${domain}.hmm - |
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
        > Ferr/${domain}.replace.tsv

    >&2 echo
done
```
