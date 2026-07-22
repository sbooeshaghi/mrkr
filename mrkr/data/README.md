# Grounding data

`gmap.txt` is the alias-inclusive human gene map used by `mrkr`. Because one alias can name more
than one gene, `gmap_canonical.txt` supplies the authoritative Ensembl gene-name pairs that take
priority during lookup. The canonical map was generated from the `gene` records in Ensembl release
114, `Homo_sapiens.GRCh38.114.gtf.gz` (GRCh38.p14; gene build updated 2024-11).

Both files are hashed into each grounded artifact. A custom two-column map remains supported and is
treated as authoritative for the organism supplied by the user.
