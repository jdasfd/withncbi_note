# HMM related resources

## PFAM-A

```bash
mkdir -p ~/data/HMM/PFAM
cd ~/data/HMM/PFAM

for basename in Pfam-A.hmm Pfam-A.hmm.dat active_site.dat; do
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/${basename}.gz
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/${basename}.gz
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/${basename}.gz
done
# wget: a non-interactive network retriever
# -N/--timestamping: don't re-retrieve files unless newer than local
# --content-disposition: honor the Content-Disposition header when choosing local file names (EXPERIMENTAL)

for basename in Pfam-A.hmm Pfam-A.hmm.dat active_site.dat; do
    echo "==> ${basename}"
    gzip -dcf ${basename}.gz > ${basename}
done
```

## TIGRFAM

```bash
mkdir -p ~/data/HMM/TIGRFAM
cd ~/data/HMM/TIGRFAM

wget -N --content-disposition ftp://ftp.jcvi.org/data/TIGRFAMs/14.0_Release/TIGRFAMs_14.0_HMM.tar.gz

mkdir -p HMM
tar --directory HMM -xzvf TIGRFAMs_14.0_HMM.tar.gz TIGR02013.HMM
tar --directory HMM -xzvf TIGRFAMs_14.0_HMM.tar.gz TIGR00485.HMM
```

## 40 single-copy genes

The protocol takes all the proteins in all the genomes under consideration and, using rapid searching and clustering algorithms, generates protein families from that complete protein set.

This systematic approach reveals 40 PhyEco marker candidates spanning the domains of bacteria and archaea.

```bash
mkdir -p ~/data/HMM/scg40
cd ~/data/HMM/scg40

curl -LO https://ndownloader.figshare.com/files/3093482

tar xvfz 3093482
```

## 120 bacterial proteins `bac120`

Ref.:
1. [Nature Microbiology, 2017, 2: 1533-1542](https://www.nature.com/articles/s41564-017-0012-7)
2. [Nature Biotechnology, 2018, 36: 996-1004](https://www.nature.com/articles/nbt.4229)

```bash
mkdir -p ~/data/HMM/bac120
cd ~/data/HMM/bac120

cp ~/Scripts/withncbi/hmm/bac120.tsv ~/data/HMM/bac120/

mkdir -p HMM

cat ~/Scripts/withncbi/hmm/bac120.tsv |
    sed '1d' |
    tsv-select -f 1 |
    grep '^TIGR' |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        tar --directory HMM -xzvf ../TIGRFAM/TIGRFAMs_14.0_HMM.tar.gz {}.HMM
    '

cat ~/Scripts/withncbi/hmm/bac120.tsv |
    sed '1d' |
    tsv-select -f 1 |
    grep -v '^TIGR' |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        curl -L http://pfam.xfam.org/family/{}/hmm > HMM/{}.HMM
    '
```
