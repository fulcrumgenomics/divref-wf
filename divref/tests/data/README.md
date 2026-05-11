# Test data

The Hail table of gnomAD variant annotations at [chr1_100001_200000.ht](chr1_100001_200000.ht) contains all 9,483 gnomAD HGDP+1KG variants in the chr1:100,001-200,000 locus from gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_variant_annotations.ht (also available at s3a://gnomad-public-us-east-1/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_variant_annotations.ht).

The Hail table of gnomAD HGDP/1KG sample metadata at [hgdp_1kg_sample_metadata.ht](hgdp_1kg_sample_metadata.ht) contains all 4,151 samples from gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.ht (also available at s3a://gnomad-public-us-east-1/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.ht), stripped to just the key (`s`) and `gnomad_population_inference` field.

The VCF file at [chr1_100001_200000.vcf.gz](chr1_100001_200000.vcf.gz) contains phased genotypes for all HGDP+1KG samples in the same locus as the Hail table above.

Because each tool depends on the output of one or more previous tools, to avoid repeating calls to the tools for each test, we generate intermediate files.

The test data was generated with a Snakemake workflow:

```bash
pixi run snakemake -j4 -s workflows/create_test_data.smk
```
