# Electron transfer related proteins

```bash
cat ASSEMBLY/Pseudomonas.assembly.pass.csv |
  tsv-summarize -H -d ',' -g RefSeq_category --count
#RefSeq_category,count
#Representative Genome,533
#,1404
#Reference Genome,15
```

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

In PIRSF database, the Pfam domain was contained, so I used the Pfam HMM.

The goal of the step is to search target protein from each strain.

### hmmsearch for Rubredoxin

```bash
cd /mnt/e/data/Pseudomonas

mkdir -p Rubr/HMM

curl -L https://pfam.xfam.org/family/PF00301/hmm > Rubr/HMM/Rubr.hmm
curl -L http://www.pantherdb.org/panther/exportHmm.jsp?acc=PTHR47627 > Rubr/HMM/PTHR47627.hmm

# using Rubredoxin.hmm to search in each strains
E_VALUE=1e-20
for domain in Rubr PTHR47627; do
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

tsv-join Rubr/Rubr.replace.tsv \
    -f Rubr/PTHR47627.replace.tsv \
    > Rubr/Rubredoxin.replace.tsv

wc -l Rubr/*.tsv
#  1997 Rubr/PTHR47627.replace.tsv
#  2089 Rubr/Rubr.replace.tsv
#  1990 Rubr/Rubredoxin.replace.tsv
# remember that 1952 strains passed were obtained

# extract all Rubr pro_seqs
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 Rubr/Rubredoxin.replace.tsv) Rubr/Rubr.fa

muscle -in Rubr/Rubr.fa -out Rubr/Rubr.aln.fa

FastTree Rubr/Rubr.aln.fa > Rubr/Rubr.aln.newick

nw_reroot Rubr/Rubr.aln.newick $(nw_labels Rubr/Rubr.aln.newick | grep -E "Bac_subti|Sta_aure") |
    nw_order -c n - \
    > Rubr/Rubr.reroot.newick

# Rubredoxin among all reference genomes
cat Rubr/Rubr.replace.tsv | grep -v 'GCF'
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
```

### Counting copy of each strain

```bash
cd /mnt/e/data/Pseudomonas

cat PROTEINS/all.strain.tsv |
    sed '1d' |
    grep -F -f <(cut -f 2 Rubr/Rubredoxin.replace.tsv) |
    cut -f 2 |
    tsv-summarize -g 1 --count \
    > Rubr/strains.copy.tsv

wc -l Rubr/strains.copy.tsv
#1505 Rubr/strains.copy.tsv
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


```

### Reference 

```bash
cd /mnt/e/data/Pseudomonas

mkdir -p Rubr/tree

Rscript -e '
library(dplyr)
library(ggtree)
gtree=read.tree("/tree/bac120.reroot.newick")
ptree=read.tree("/Rubr/Rubr.reroot.newick")
p1 <- ggtree(gtree)
p2 <- ggtree(ptree)

d1 <- p1$data
d2 <- p2$data

d2$x
'

library(dplyr)
library(ggtree)
tree1=read.tree("bac120.model.reroot.newick")
tree2=read.tree("bac120.protein.newick")
p1 <- ggtree(tree1)
p2 <- ggtree(tree2)

d1 <- p1$data
d2 <- p2$data

d2$x <- max(d2$x) - d2$x + max(d1$x) + 1 #翻转第二棵树
d2$y <- d2$y + 8 #将第二棵树向上移动对齐
p3 <- ggtree::rotate(p1,29) #旋转node，使两棵树拓扑结构一致

pp <- p3 + geom_tiplab(offset=0.05) + geom_treescale()+ geom_highlight(node=25,fill="red")+ geom_tree(data=d2) + geom_tiplab(data = d2, hjust=1, offset =-0.05)+geom_treescale(x=2)

#geom_tiplab(offset=0.05) 添加label，offset：label与树枝末端距离
#geom_treescale() #添加表尺
#geom_highlight(node=25,fill="red") #添加红色区域
```

## Compare protein seqs with TIGRFAM

TIGRFAM database could be seen in [hmm.md](hmm.md)

```bash
cd /mnt/e/data/Pseudomonas

cat ~/data/HMM/TIGRFAM/HMM/*.HMM > Rubr/HMM/TIGRFAM.hmm
cat ~/data/HMM/PGAP/HMM/hmm_PGAP/*.HMM > Rubr/HMM/PGAP.hmm

# hmmpress to prepare an database
hmmpress Rubr/HMM/TIGRFAM.hmm
hmmpress Rubr/HMM/PGAP.hmm
# hmmpress could prepare an HMM database for faster hmmscan searches

# faops extract protein seq
faops some PROTEINS/all.replace.fa <(tsv-select -f 2 Rubr/Rubredoxin.replace.tsv) Rubr/Rubr.fa

# using protein seqs to scan directly in TIGRFAM database
E_VALUE=1e-10
for NAME in Rubr; do
    >&2 echo "==> Protein [${NAME}]"
        if [ ! -s ${NAME}/HMM/TIGRFAM.hmm ]; then
        echo no TIGRFAM
        exit
        fi
    
    hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw \
        -o ${NAME}/${NAME}.tigrfam.txt --tblout ${NAME}/${NAME}.tigrfam.tbl \
        ${NAME}/HMM/TIGRFAM.hmm ${NAME}/${NAME}.fa
    
    echo hmmscan complete
done

# using protein seqs to scan directly in PGAP database (included TIGR)
E_VALUE=1e-10
for NAME in Rubr; do
    >&2 echo "==> Protein [${NAME}]"
        if [ ! -s ${NAME}/HMM/PGAP.hmm ]; then
        echo no PGAP
        exit
        fi
    
    hmmscan --cpu 6 -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw \
        -o ${NAME}/${NAME}.pgap.txt --tblout ${NAME}/${NAME}.pgap.tbl \
        ${NAME}/HMM/PGAP.hmm ${NAME}/${NAME}.fa
    
    echo hmmscan complete
done

# reformat and extract results from tbl
cat Rubr/Rubr.tigrfam.tbl |
    grep '^TIGR' |
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
    ' > Rubr/Rubr.tigrfam.tsv

cat Rubr/Rubr.tigrfam.tsv |
    tsv-filter --le 4:1e-50 |
    wc -l
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
