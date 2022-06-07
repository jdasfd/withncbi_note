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

## IPR000813 - 7Fe ferredoxin

- PRINTS: PR00354 7Fe ferredoxin signature

PRINTS database does not have the HMM of the protein.

7FE8SFRDOXIN is a 3-element fingerprint that provides a signature for 7Fe ferredoxins. The fingerprint was derived from an initial alignment of 2 sequences: motif 1 contains 2 of 3 invariant Cys residues contributing to the 3Fe-4S cluster, which is located upstream from the 4Fe-4S domain; motif 2 contains the first invariant Cys of the 4Fe-4S cluster; and motif 3 contains 3 invariant cysteines of the 4Fe-4S cluster and one Cys ligand of its downstream 3Fe-4S counterpart (see scheme above).

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
```

## IPR024922 - Rubredoxin

- PIRSF: PIRSF000071 Rubredoxin
- Pfam: PF00301: Rubredoxin

In PIRSF database, the Pfam domain was contained, so I used the Pfam HMM.

```bash
cd /mnt/e/data/Pseudomonas

mkdir -p Rubr/HMM

curl -L https://pfam.xfam.org/family/PF00301/hmm > Rubr/HMM/Rubredoxin.hmm

E_VALUE=1e-20
for domain in Rubredoxin; do
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
```
