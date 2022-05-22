# Aligning various genera missing from `bacteria_gr.md`

The markdown was originated from my teacher [bacteria_gr.md](https://github.com/wang-q/withncbi/blob/master/taxon/bacteria_gr.md). Noted and recorded my own understanding of those command.

## Strain info

```bash
export RANK_NAME=bacteria_ar

mkdir -p /mnt/e/data/alignment/${RANK_NAME}
cd /mnt/e/data/alignment/${RANK_NAME}


```

```bash
echo '
.headers ON

SELECT
    organism_name,
    species,
    genus,
    ftp_path,
    assembly_level
FROM ar
WHERE 1=1
    AND tax_id != species_id               -- with strain ID
    AND species_id IN (29388)
' |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > Scap.assembly.tsv

echo '
SELECT
    species || " " || REPLACE(assembly_accession, ".", "_") AS organism_name,
    species,
    genus,
    ftp_path,
    assembly_level
FROM ar
WHERE 1=1
    AND tax_id = species_id               -- no strain ID
    AND assembly_level IN ("Chromosome", "Complete Genome")
    AND species_id IN (29388)
' |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >> Scap.assembly.tsv

echo '
SELECT
    species || " " || REPLACE(assembly_accession, ".", "_") AS organism_name,
    species,
    genus,
    ftp_path,
    assembly_level
FROM ar
WHERE 1=1
    AND tax_id = species_id               -- no strain ID
    AND assembly_level IN ("Chromosome", "Complete Genome")
' |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > SPECIES_ID.txt
```
