# Learning and Classifying Plasmids

My script was originated from my teacher [plasmid.md](https://github.com/wang-q/withncbi/blob/master/taxon/plasmid.md).

I rewrote some commands according to my own demanding and added some notes in correspoding locations.

## NCBI RefSeq

```bash
mkdir -p /mnt/d/data/plasmid
cd /mnt/d/data/plasmid

rsync -avP ftp.ncbi.nlm.nih.gov::refseq/release/plasmid/ RefSeq/
```

- Get taxon id and locus name.

```bash
gzip -dcf RefSeq/*.genomic.gbff.gz > genomic.gbff
# -d, --decompress: decompress
# -c, --stdout: write on standard output, keep original files unchanged
# -f, --force: force overwrite of output file and compress links
# so this step could decompress all plasmids' gbff.gz into a genomic.gbff file

perl ~/Scripts/withncbi/taxon/gb_taxon_locus.pl genomic.gbff > refseq_id_seq.csv
#There are [43857] sequences.
#There are [43857] valid sequences.
rm genomic.gbff
```

- Retrive fasta and check.

```bash
gzip -dcf RefSeq/plasmid.1.1.genomic.fna.gz |
    grep "^>" |
    head -n 5
#>NZ_D13972.1 Synechococcus sp. PCC 7002 strain PR-6 plasmid Plasmid pAQ1, complete sequence
#>NZ_Y18549.1 Escherichia coli strain K5533 plasmid pColK-JA533, complete sequence
#>NZ_U51470.1 Pasteurella multocida plasmid unnamed, complete sequence
#>NZ_M10917.1 Bacillus thuringiensis plasmid unnamed, complete sequence
#>NZ_U40997.1 Listeria monocytogenes strain BM4293 plasmid pIP823, complete sequence
```

So we got the plasmid taxon in `>NZ_...`, it the same with `refseq_id_seq.csv` col1.

- faops counting N50

The N metrics are a measure of contiguity of a set of sequences often used to assess genome assemblies. The N50 is related to the median and mean length of a set of sequences. [Ref](https://timkahlke.github.io/LongRead_tutorials/APP_MET.html)

```bash
faops n50 -S -C RefSeq/*.genomic.fna.gz
#N50     209488
#S       4158684992
#C       43857

# -S: compute sum of size of all entries
# -C: count entries

gzip -dcf RefSeq/*.genomic.fna.gz > RefSeq/plasmid.fa
```

## MinHash to get non-redundant plasmids

```bash
mkdir /mnt/d/data/plasmid/nr
cd /mnt/d/data/plasmid/nr

faops size ../RefSeq/plasmid.fa > refseq.sizes

tsv-filter refseq.sizes --le 2:2000 | wc -l
#3180

faops some ../RefSeq/plasmid.fa <(tsv-filter refseq.sizes --gt 2:2000) refseq.fa
# extract plasmids over 2000 from plasmid.fa to refseq.fa
```

- MinHash algorithm for representing reads

[MinHash](https://en.wikipedia.org/wiki/MinHash) (or the min-wise independent permutations locality sensitive hashing scheme) is a technique for quickly estimating how similar two sets are.

Sketches could reduce representations of sequences for `mash` comparison. Sketches allow fast distance estimations with low storage and memory requirements.

[Explanation](https://mash.readthedocs.io/en/latest/sketches.html): Each k-mer in a sequence is hashed, which creates a pseudo-random identifier. By sorting them, a small subset can represent the entire sequence (min-hashes). The more similar another sequence is, the more min-hashes it is likely to share.

`mash sketch -h` could give you the help information.

```bash
cd /mnt/d/data/plasmid/nr

cat refseq.fa | 
    mash sketch -k 21 -s 1000 -i -p 8 - -o refseq.plasmid.k21s1000.msh
# -k: k-mer size
# -s: sketch size, correponds to the number of (non-redundant) min-hashes that are kept
# -i: sketch individual sequences, rather than whole files
# -p: threads
# -: read from standard input
```

- Split

`split`: output pieces of FILE to PREFIXaa, PREFIXab, ...; default size is 1000 lines, and default PREFIX is 'x'.

`split [OPTION]... [FILE [PREFIX]]`

`find`: find files under the dir

`find [-H] [-L] [-P] [-Olevel] [-D debugopts] [path...] [expression]`

This step is using `split` to split a big file `refseq.fa` into small pieces, then using files splited in `mash sketch` (parallel could reduce time).

```bash
cd /mnt/d/data/plasmid/nr
mkdir -p job
faops size refseq.fa |
    cut -f 1 |
    split -l 1000 -a 3 -d - job/
# split arguments:
# -l, --lines=NUMBER: generate CHUNKS output files
# CHUNKS: N - split into N files based on size of input;
#         K/N - output Kth of N to stdout
#         l/N - split into N files without splitting lines/records
#         l/K/N - output Kth of N to stdout without splitting lines/records
#         r/N - like 'l' but use round robin distribution
#         r/K/N - likewise but only output Kth of N to stdout
# -a: generate suffixes of length N (default 2)
# -d: use numeric suffixes starting at 0, not alphabetic
# -a 3 -d means the output file will be: 000, 001, 002, 003
# -: read from stdin

cat job/000 | wc -l
#1000
# The refseq.fa were splited into multiple files with 1000 taxons a file

# MinHash for all files sketch
find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        faops some refseq.fa {} stdout |
            mash sketch -k 21 -s 1000 -i -p 6 - -o {}.msh
    '
# -type: f: a regular file; d: directory; l: symbolic link; c: character devices; b: block devices; p: named pipe; s: socket
# -maxdepth [a non-negative integer]: Descend at most levels of directories below the command line arguments.
# --line-buffer: buffer output on line basis.
```

- Dist among all .msh sketch files

Estimate the distance of each query sequence to the reference. Both the reference and queries can be Mash sketch files (`.msh`) with matching k-mer sizes. The output fields are: (reference-ID, query-ID, distance, p-value, shared-hashes).

`mash dist -h` could give you the help information.

Usage: `mash dist [options] <reference> <query> [<query>] ...`

```bash
cd /mnt/d/data/plasmid/nr

# Count distance from all .msh results
find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        mash dist -p 6 {}.msh refseq.plasmid.k21s1000.msh > {}.tsv
    '
# -p <int>: threads

cat job/000.tsv | head
#NZ_D13972.1     NZ_D13972.1     0       0       1000/1000
#NZ_M10917.1     NZ_D13972.1     1       1       0/1000
#NZ_U40997.1     NZ_D13972.1     1       1       0/1000
#NZ_U35036.1     NZ_D13972.1     1       1       0/1000
#NZ_L25424.1     NZ_D13972.1     1       1       0/1000
#NZ_M60875.1     NZ_D13972.1     1       1       0/1000
#NC_021737.1     NZ_D13972.1     1       1       0/1000
#NC_010097.1     NZ_D13972.1     1       1       0/1000
#NC_002100.1     NZ_D13972.1     1       1       0/1000
#NZ_CP010448.1   NZ_D13972.1     1       1       0/1000
#each col means:
#ref_ID query_ID    dist    p-val   shared_hashes
```

- Use dist to find all redundant plasmids

```bash
# distance < 0.01
find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 16 '
        cat {}.tsv |
            tsv-filter --ff-str-ne 1:2 --le 3:0.01
    ' \
    > redundant.tsv
# tsv-filter:
# --ff-str-ne 1:2 : means strings in col1 and col2 not equal
# --le 3:0.01 : means col3 smaller or equal 0.01
# from the above, which means distance bewtween ref_seq and query_seq was smaller than 0.01

head -n 5 redundant.tsv
#NZ_KX777254.1   NC_003277.2     0.000167546     0       993/1000
#NZ_KU761328.1   NZ_CP011983.1   0.00570819      0       797/1000
#NZ_CP080361.1   NZ_CP011983.1   0.00550978      0       803/1000
#NZ_CP031725.1   NZ_CP011983.1   0.00567497      0       798/1000
#NZ_KX777254.1   NZ_CP039586.1   0.000215742     0       991/1000

cat redundant.tsv | wc -l
#950476
```

- Deal with those redundant plasmids

```bash
cat redundant.tsv |
    perl -nla -F"\t" -MGraph::Undirected -e '
        BEGIN {
            our $g = Graph::Undirected->new;
        }

        $g->add_edge($F[0], $F[1]);

        END {
            for my $cc ( $g->connected_components ) {
                print join qq{\t}, sort @{$cc};
            }
        }
    ' \
    > connected_components.tsv

# Perl One-liner:
# -n: assume while (<>) { ... } loop around program
# -p: assume loop like -n but print line also, like sed
# -a: autosplit mode with -n or -p (splits $_ into @F)
# -F: split() pattern for -a switch (//'s are optional)
# -l[octal]: enable line ending processing, specifies line terminator
# -[mM][-]module: execute use/no module... before executing program
# Graph::Undirected allows you to create undirected graphs.
```
