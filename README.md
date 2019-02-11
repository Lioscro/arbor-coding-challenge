# arbor-coding-challenge
Repository for the Arbor Biotechnologies internship coding challenge. This
README describes the requirements, usage, and output of
[extract_protein_info.py](src/extract_protein_info.py). In a general sense,
this command-line utility reads in a GenBank format .gbff file, parses its
CDS information, and outputs useful information about the proteins (CDS) in
the file. It is also able to fetch predicted protein conserved domains by
querying the NCBI Conserved Domains database via POST request.

## Requirements
The versions indicate the versions of the packages that were installed on
the machine on which the script was tested on.
- Python 3.7.1
- [Requests](https://github.com/kennethreitz/requests) 2.21.0
- [Numpy](http://www.numpy.org/) 1.15.4
- [Pandas](https://pandas.pydata.org/) 0.23.4
- [Biopython](https://biopython.org/) 1.72

## Usage
```
usage: extract_protein_info.py [-h] [-p P] [-g GBFF_FILTER] [-d DOMAIN_FILTER]
                               [-r] [-o]
                               gbff

positional arguments:
  gbff                  .gbff file containing the protein(s) of interest

optional arguments:
  -h, --help            show this help message and exit
  -p P                  prefix for the output file(s)
  -g GBFF_FILTER, --gbff-filter GBFF_FILTER
                        a comma-separated list of proteins to filter for; can
                        be protein_id, locus_tag, translation, product
  -d DOMAIN_FILTER, --domain-filter DOMAIN_FILTER
                        a comma-separated list of strings to filter proteins
                        by their conserved domain; applied only when domains
                        are retrieved
  -r, --repeat-only     only consider proteins that have < infinity distance
                        from a repeat_region
  -o, --override        use this option to always fetch conserved protein
                        domains--MAY TAKE A LONG TIME
```

## Features
This script has a few notable features.
- Calculates the distance to nearby repeat regions. If no repeat region is
found in the same record, then it is given a distance of `Numpy.Infinity`.
- Submits queries to the
[NCBI Conserved Domain database](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml)
to retrieve predicted conserved domains in the proteins. This feature is best
used with either `--gbff-filter` or `--repeat-only` options. Otherwise, the
domain prediction may take a long time.
- A two-step filtering process. If the `--gbff-filter` option is used, the
script first filters the input GenBank .gbff file by the `protein_id`,
`locus_tag`, `translation` or `product` of each CDS. If the `--repeat-only`
flag is used, only proteins with a nearby repeat_region will be considered
(i.e. proteins with a repeat_region in the same record). If the `--domain-filter`
is also given, assuming that the domains were successfully retrieved, filters
for proteins (CDS) that have all the filters in either the `Superfamily` or
`Short name` fields.

## Output
Depending on whether or not filters were applied or the `--override` flag was
used, this script produces one or two files. The script always outputs
a table of filtered (or non-filtered) proteins in CSV format. If protein
domains were retrieved, then it also outputs a second CSV file with
domain predictions of proteins with `locus_tag` corresponding to those of the
first output.

## Filtering for two instances of Cas9
The robust filtering process implemented by this script makes it relatively
simply to find proteins with domains of interest. The two instances of Cas9,
corresponding to protein ids `WP_053019794.1` and `WP_010922251.1` in the two
.gbff files in the [sample_gbff](sample_gbff) folder can be found with the
following commands.
`python extract_protein_info.py GCF_001239625.1_7068_7_24_genomic.gbff --repeat-only -d Cas9`
and
`python extract_protein_info.py GCF_002014815.1_ASM201481v1_genomic.gbff --repeat-only -d Cas9`
In both cases, the script first filters the gbff file for proteins
that are close (distance < `Numpy.Infinity`) to repeat regions,
retrieves conserved domain predictions of these proteins, and then
filters the table of conserved domains for the string `Cas9`.
