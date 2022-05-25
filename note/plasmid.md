# Learning and Classifying Plasmids

My script was originated from my teacher [plasmid.md](https://github.com/wang-q/withncbi/blob/master/taxon/plasmid.md).

I rewrote some commands according to my own demanding and added some notes in correspoding locations.

- [Learning and Classifying Plasmids](#learning-and-classifying-plasmids)
  - [NCBI RefSeq](#ncbi-refseq)
  - [MinHash to get non-redundant plasmids](#minhash-to-get-non-redundant-plasmids)
  - [Grouping by MinHash](#grouping-by-minhash)
  - [Plasmid: prepare](#plasmid-prepare)
  - [Plamid: run](#plamid-run)

## NCBI RefSeq

```bash
mkdir -p /mnt/d/data/plasmid
cd /mnt/d/data/plasmid

rsync -avP ftp.ncbi.nlm.nih.gov::refseq/release/plasmid/RefSeq/
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

My understanding of this:

Because of the redundant.tsv contained info of all nodes distance that smaller than 0.01.

So using the `Graph::Undirected` to link all nodes. Those repeated nodes will be excluded.

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
# Graph::Undirected:
# allows you to create undirected graphs
# add_edge: Add the edge to the graph
# connected_components: for an undirected graph, returns the vertices of the connected components of the graph as a list of anonymous arrays

cat connected_components.tsv | head -n 3
#NZ_CP059740.1   NZ_CP059748.1   NZ_CP059794.1   NZ_CP059817.1   NZ_CP062266.1   NZ_CP072708.1
#NC_016979.1     NZ_CP023921.1   NZ_CP024043.1   NZ_CP029589.1   NZ_CP029592.1   NZ_CP045021.1   NZ_CP054976.1  NZ_CP069457.1   NZ_CP069492.1   NZ_CP069555.1   NZ_CP069585.1   NZ_CP079119.1   NZ_CP079177.1   NZ_CP079622.1   NZ_CP079650.1  NZ_CP079662.1    NZ_CP079696.1   NZ_CP079812.1   NZ_CP083065.1   NZ_LR890626.1
#NZ_CP021666.1   NZ_CP022935.1

cat connected_components.tsv |
    perl -nla -F"\t" -e 'printf qq{%s\n}, $_ for @F' \
    > components.list
# change all components to one name a line

wc -l connected_components.tsv components.list
#3920 connected_components.tsv
#24159 components.list
#28079 total
```

[Graph (discrete mathematics)](https://en.wikipedia.org/wiki/Graph_(discrete_mathematics)):

Undirected graph:  if the vertices represent people at a party, and there is an edge between two people if they shake hands, then this graph is undirected because any person A can shake hands with a person B only if B also shakes hands with A.

Directed graph: if any edge from a person A to a person B corresponds to A owes money to B, then this graph is directed, because owing money is not necessarily reciprocated.

- Extract those redundant plasmids

```bash
faops some -i refseq.fa components.list stdout > refseq.nr.fa
faops some refseq.fa <(cut -f 1 connected_components.tsv) stdout >> refseq.nr.fa

rm -fr job
```

## Grouping by MinHash

All non-redundant plasmids got from the previous steps will be grouped by MinHash again.

- MinHash algorithm for representing reads

```bash
mkdir /mnt/d/data/plasmid/grouping
cd /mnt/d/data/plasmid/grouping

cat ../nr/refseq.nr.fa |
    mash sketch -k 21 -s 1000 -i -p 8 - -o refseq.nr.k21s1000.msh
```

- Split

```bash
mkdir -p job
faops size ../nr/refseq.nr.fa |
    cut -f 1 |
    split -l 1000 -a 3 -d - job/

find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        faops some ../nr/refseq.nr.fa {} stdout |
            mash sketch -k 21 -s 1000 -i -p 6 - -o {}.msh
    '

find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        mash dist -p 6 {}.msh refseq.nr.k21s1000.msh > {}.tsv
    '

find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 1 '
        cat {}.tsv
    ' \
    > dist_full.tsv
```

- Use dist to find all connected plasmids

```bash
# distance < 0.05
cat dist_full.tsv |
    tsv-filter --ff-str-ne 1:2 --le 3:0.05 \
    > connected.tsv

head -n 5 connected.tsv
#NZ_CP035352.1   NG_048225.1     0.047688        0       225/1000
#NZ_CP025252.1   NZ_CP011982.1   0.0415686       0       264/1000
#NZ_LT985257.1   NZ_CP011982.1   0.0275126       0       390/1000
#NZ_CP053253.1   NZ_CP014966.1   0.0471738       0       228/1000
#NZ_CP018774.2   NZ_CP014966.1   0.0458408       0       236/1000

cat connected.tsv | wc -l
#179616
```

- Group according to dist

```bash
mkdir -p group

cat connected.tsv |
    perl -nla -F"\t" -MGraph::Undirected -MPath::Tiny -e '
        BEGIN {
            our $g = Graph::Undirected->new;
        }

        $g->add_edge($F[0], $F[1]);

        END {
            my @rare;
            my $serial = 1;
            my @ccs = $g->connected_components;
            @ccs = map { $_->[0] }
                sort { $b->[1] <=> $a->[1] }
                map { [ $_, scalar( @{$_} ) ] } @ccs;
            for my $cc ( @ccs ) {
                my $count = scalar @{$cc};
                if ($count < 50) {
                    push @rare, @{$cc};
                }
                else {
                    path(qq{group/$serial.lst})->spew(map {qq{$_\n}} @{$cc});
                    $serial++;
                }
            }
            path(qq{group/00.lst})->spew(map {qq{$_\n}} @rare);

            path(qq{grouped.lst})->spew(map {qq{$_\n}} $g->vertices);
        }
    '
```

> **Explain code**:
>
> - Schwartzian transform:
>
> ```perl
> @sorted = map  { $_->[0] }
>          sort { $a->[1] <=> $b->[1] or $a->[0] cmp $b->[0] } # Use numeric comparison, fall back to string sort on original
>          map  { [$_, length($_)] }    # Calculate the length of the string
>               @unsorted;
> ```
>
> - Path::tiny
>
> ```perl
> use Path::Tiny;
> 
> # creating Path::Tiny objects
> 
> $dir = path("/tmp");
> $foo = path("foo.txt");
>
> $subdir = $dir->child("foo");
> $bar = $subdir->child("bar.txt");
>
> # writing files
> 
> $bar->spew( @data );
> $bar->spew_utf8( @data );
> ```
>
> So the code in:
>
> ```perl
> path(qq{group/00.lst})->spew(map {qq{$_\n}} @rare);
> ```
>
> means write @rare each element a line into the ./group/00.lst file.

- Get non-grouped plasmids

```bash
# get non-grouped
# this will no be divided to subgroups
faops some -i ../nr/refseq.nr.fa grouped.lst stdout |
    faops size stdin |
    cut -f 1 \
    > group/lonely.lst

wc -l group/*
#  5311 group/00.lst
#  4208 group/1.lst
#   695 group/2.lst
#   189 group/3.lst
#    97 group/4.lst
#    97 group/5.lst
#    86 group/6.lst
#    64 group/7.lst
#    64 group/8.lst
#    54 group/9.lst
#  9583 group/lonely.lst
# 20448 total
```

- Analyze each group by R

```bash
find group -maxdepth 1 -type f -name "[0-9]*.lst" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"

        faops some ../nr/refseq.nr.fa {} stdout |
            mash sketch -k 21 -s 1000 -i -p 6 - -o {}.msh

        mash dist -p 6 {}.msh {}.msh > {}.tsv
    '

find group -maxdepth 1 -type f -name "[0-9]*.lst.tsv" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"

        cat {} |
            tsv-select -f 1-3 |
            Rscript -e '\''
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
                write.tree(phy=tree, file="{.}.tree.nwk")

                group <- cutree(clusters, h=0.2) # k=3
                groups <- as.data.frame(group)
                groups$ids <- rownames(groups)
                rownames(groups) <- NULL
                groups <- groups[order(groups$group), ]
                write_tsv(groups, "{.}.groups.tsv")
            '\''
    '
```

> **Explain code**:
>
> `read_tsv` will give the below info:
>
> ```txt
> -- Column specification >-----------------------------------------------
> Delimiter: "\t"
> chr (2): X1, X2
> dbl (1): X3
> ```
>
> `ape` means Analysis of Phylogenetics and Evolution R packages.
>
> `pivot_wider()` "widens" data, increasing the number of columns and decreasing the number of rows. This step will change the `pair_dist` to a matrix - each col is from `pair_dist$X2` and values are from `pair_dist$X3`.
>
> The main goal is to transform all pair values into a matrix.
>
> `tmp[,-1]` means accessing tmp all cols without the 1st name col.
>
> `tmp[,1]` means accessing tmp 1st col as rownames.
>
> `ward.D2` - Like most other clustering methods, Ward’s method is computationally intensive. However, Ward’s has significantly fewer computations than other methods.
>
> `hclust` - [what is hierarchical clustering](https://www.displayr.com/what-is-hierarchical-clustering/)
>
> `cutree` - Cut a Tree into Groups of Data. Cuts a tree, *e.g.*, as resulting from hclust, into several groups either by specifying the desired number(s) of groups or the cut height(s).
>
> After all those steps in R, we converted pair dist tsv files into dist matrices, which were used for each group building a phylogenic tree and were cutted into groups based on the tree.

- Analyze subgroup

```bash
# subgroup
mkdir -p subgroup
cp group/lonely.lst subgroup/

find group -name "*.groups.tsv" | sort |
    parallel -j 1 -k '
        cat {} | sed -e "1d" | xargs -I[] echo "{/.}_[]"
    ' |
    sed -e 's/.lst.groups_/_/' |
    perl -na -F"\t" -MPath::Tiny -e '
        path(qq{subgroup/$F[0].lst})->append(qq{$F[1]});
    '
# -k/--keep-order: keep sequence of output same as the order of input
# {/.}: Basename of input line without extension

# ignore small subgroups
find subgroup -name "*.lst" | sort |
    parallel -j 1 -k '
        lines=$(cat {} | wc -l)

        if (( lines < 5 )); then
            echo -e "{}\t$lines"
            cat {} >> subgroup/lonely.lst
            rm {}
        fi
    '
# (( ... )) this performs arithmetic
# Like [[ ... ]], it returns an exit code of zero (true) if the result nonzero

# append ccs
# So this step is actually appending connected_components to a subgroup file
cat ../nr/connected_components.tsv |
    parallel -j 1 --colsep "\t" '
        file=$(rg -F -l  "{1}" subgroup)
        echo {} | tr "[:blank:]" "\n" >> ${file}
    '
# rg: recursively search the current directory for lines matching a pattern
# -F/--fixed-strings: treat the pattern as a literal string instead of a regular expression
# -l/--file-with-matches: print the paths with at least one match and suppress match contents
# tr: copies the standard input to the standard output with substitution or deletion of selected characters
# "[:blank:]": all horizontal whitespace - POSIX Basic Regular Expressions
```

> **Explanation**:
>
> `parallel {}`:
>
> `{}` means input line. *E.g.*, `group/1.lst.groups.tsv` was passed to parallel from stdout, so the `cat {}` actually meant `cat group/1.lst.groups.tsv`.
>
> `{.}` means input line without extension. So `{.}` actually meant `group/1.lst.groups` (without `.tsv`).
>
> `{/}` means basename of input line removed. So `{/}` actually meant `1.lst.groups.tsv` (without `group/`).
>
> `{//}` means dirname. So `{//}` actually meant `group` (without `/1.lst.groups.tsv`).
>
> `{/.}` means basename of input line without extension. So `{/.}` actually meant `1.lst.groups`.
>
> `{#}` means sequence number of the job to run. So `{#}` actually meant `1`.
>
> `parallel {1}`: will use an example for explaining.
>
> ```bash
> parallel echo {1} {2} {3} ::: 6 7 ::: 4 5 ::: 1 2 3 | head -n 5
> #6 4 1
> #6 4 2
> #6 4 3
> #6 5 1
> #6 5 2
> # multiple input will give out an combination input like (6,4,1), (6,4,2), (6,4,3) ... (7,5,2), (7,5,3)
> # then {1} means col1, {2} means col2
>
> parallel echo {2} ::: 6 7 ::: 4 5 ::: 1 2 3 | head -n 5
> #4
> #4
> #4
> #5
> #5
> # only output col2
> ```
>
> `--colseq/-C regexp`: column separator
>
> ```bash
> parallel echo {4} {3} {2} {1} \
> ::: A-B C-D ::: e-f g-h
> #e-f A-B
> #g-h A-B
> #e-f C-D
> #g-h C-D
> # without col4 and col3, so the result was the same to the below one
>
> parallel echo {2} {1} ::: A-B C-D ::: e-f g-h
> #e-f A-B
> #g-h A-B
> #e-f C-D
> #g-h C-D
>
> parallel --colsep '-' echo {4} {3} {2} {1} \
> ::: A-B C-D ::: e-f g-h
> #f e B A
> #h g B A
> #f e D C
> #h g D C
> # A B was seperated by -, so the table would be (A B e f) for the row1, (A B g h) for the row2
>
> parallel --colsep '-' echo {2} {1} \
> ::: A-B C-D ::: e-f g-h
> #B A
> #B A
> #D C
> #D C
> ```

```bash
# remove duplicates
find subgroup -name "*.lst" | sort |
    parallel -j 1 '
        cat {} | sort | uniq > tmp.lst
        mv tmp.lst {}
    '

wc -l subgroup/* |
    sort -nr |
    head -n 100

wc -l subgroup/* |
    perl -pe 's/^\s+//' |
    tsv-filter -d" " --le 1:10 |
    wc -l
#251

wc -l subgroup/* |
    perl -pe 's/^\s+//' |
    tsv-filter -d" " --ge 1:50 |
    tsv-filter -d " " --regex '2:\d+' |
    sort -nr \
    > next.tsv

wc -l next.tsv
#104

# rm -rf job
```

## Plasmid: prepare

- Split sequences

```bash
mkdir /mnt/d/data/plasmid/GENOMES
mkdir /mnt/d/data/plasmid/taxon

cd /mnt/d/data/plasmid/grouping

echo -e "#Serial\tGroup\tCount\tTarget" > ../taxon/group_target.tsv

cat next.tsv |
    cut -d" " -f 2 |
    parallel -j 4 -k --line-buffer '
        echo >&2 "==> {}"

        GROUP_NAME={/.}
        TARGET_NAME=$(head -n 1 {} | perl -pe "s/\.\d+//g")

        SERIAL={#}
        COUNT=$(cat {} | wc -l)

        echo -e "${SERIAL}\t${GROUP_NAME}\t${COUNT}\t${TARGET_NAME}" >> ../taxon/group_target.tsv

        faops order ../nr/refseq.fa {} stdout |
            faops filter -s stdin stdout \
            > ../GENOMES/${GROUP_NAME}.fa
    '

cat next.tsv |
    cut -d" " -f 2 |
    parallel -j 4 -k --line-buffer '
        echo >&2 "==> {}"
        GROUP_NAME={/.}
        faops size ../GENOMES/${GROUP_NAME}.fa > ../taxon/${GROUP_NAME}.sizes
    '

# split-name
find ../GENOMES -maxdepth 1 -type f -name "*.fa" | sort |
    parallel -j 4 '
        faops split-name {} {.}
    '
# faops split-name: Split an fa file into several files/Using sequence names as file names

# mv to dir of basename
find ../GENOMES -maxdepth 2 -mindepth 2 -type f -name "*.fa" | sort |
    parallel -j 4 '
        mkdir -p {.}
        mv {} {.}
    '
# because GENOMES had *.fa and dir with splited *.fa, so the arguments should be specified
# that is -maxdepth 2 -mindepth 2: which means find the GENOMES subdir *.fa
```

> **Explanation**:
>
> With `--keep-order --line-buffer` will output lines from the first job continuously while it is running, then lines from the second job while that is running. It will buffer full lines, but jobs will not mix.
>
> `next.tsv` gives out the `strains_num subgroup/<file>.lst`, which means different subgroup and its strains. Using the `cut` to give the file name into parallel, each line a file. `TARGET_NAME` was the first strain in `*.lst`. In `*.lst` file, all strains accession were contained for `faops filter` extraction its genome info into dir GENOMES.
>
> The split-name step was using `faops split-name [options] <in.fa> <outdir>`. The genomes of each plasmid would be seperated into a dir named as `*.fa`/`*.lst` but without extensions.
>
> Then move the subdir `*.fa` into its own name dir.

- `prepseq`

```bash
cd /mnt/d/data/plasmid/

cat taxon/group_target.tsv |
    sed -e '1d' |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
        echo -e "==> Group: [{2}]\tTarget: [{4}]\n"

        for name in $(cat taxon/{2}.sizes | cut -f 1); do
            egaz prepseq GENOMES/{2}/${name}
        done
    '
# e.g. a line in group_target.tsv | sed -e '1d' : "1       1_4     679     NC_002122"
# so that the next command line could be translated in:
# for name in $(cat taxon/1_4.sizes | cut -f 1); do
#   egaz prepseq GENOMES/1_4/NC_002122 ...
# done
# because we put the strains.fa into its own name dir, and egaz prepseq will accept the dir
```

- Check outliers of lengths

```bash
cd /mnt/d/data/plasmid/

cat taxon/*.sizes | cut -f 1 | wc -l
#14454

cat taxon/*.sizes | cut -f 2 | paste -sd+ | bc
#1584814437
# paste:
# -s (serial): reads all the lines from a single file and merges all these lines into a single line with each line separated by tab
# -d (delimiter): tab as default, this can change any other character
# bc: used for command line calculator

cat taxon/group_target.tsv |
    sed -e '1d' |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
        echo -e "==> Group: [{2}]\tTarget: [{4}]"

        median=$(cat taxon/{2}.sizes | datamash median 2)
        mad=$(cat taxon/{2}.sizes | datamash mad 2)
        lower_limit=$( bc <<< " (${median} - 2 * ${mad}) / 2" )

#        echo $median $mad $lower_limit
        lines=$(tsv-filter taxon/{2}.sizes --le "2:${lower_limit}" | wc -l)

        if (( lines > 0 )); then
            echo >&2 "    $lines lines to be filtered"
            tsv-join taxon/{2}.sizes -e -f <(
                    tsv-filter taxon/{2}.sizes --le "2:${lower_limit}"
                ) \
                > taxon/{2}.filtered.sizes
            mv taxon/{2}.filtered.sizes taxon/{2}.sizes
        fi
    '
# tsv-join -e: exclude matching records

cat taxon/*.sizes | cut -f 1 | wc -l
#14414

cat taxon/*.sizes | cut -f 2 | paste -sd+ | bc
#1583506331
```

- Rsync to HPCC

```bash
rsync -avP /mnt/d/data/plasmid \
    wangq@202.119.37.251:jyq/data/
# if /mnt/d/data/plasmid do not contained the '/' at the end, then it will transform a whole dir to hpcc\
```

## Plamid: run

After sending the dir to HPCC, the operating path will be the path on HPCC.

```bash
cd jyq/data/plasmid/

cat taxon/group_target.tsv |
    sed -e '1d' | grep "^53" |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo -e "==> Group: [{2}]\tTarget: [{4}]\n"

        egaz template \
            GENOMES/{2}/{4} \
            $(cat taxon/{2}.sizes | cut -f 1 | grep -v -x "{4}" | xargs -I[] echo "GENOMES/{2}/[]") \
            --multi -o groups/{2}/ \
            --order \
            --parallel 24 -v

#        bash groups/{2}/1_pair.sh
#        bash groups/{2}/3_multi.sh

        bsub -q mpi -n 24 -J "{2}-1_pair" "bash groups/{2}/1_pair.sh"
        bsub -w "ended({2}-1_pair)" \
            -q mpi -n 24 -J "{2}-3_multi" "bash groups/{2}/3_multi.sh"
    '

# clean
find groups -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find groups -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type f -name "output.*" | parallel -r rm

echo \
    $(find groups -mindepth 1 -maxdepth 1 -type d | wc -l) \
    $(find groups -mindepth 1 -maxdepth 3 -type f -name "*.nwk.pdf" | wc -l)
```
