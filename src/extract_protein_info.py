'''
extract_protein_info.py

This script is a command-line utility that parses an input .gbff file
and parses the proteins (CDS) in the file according to any filters or
command-line arguments that were provided (i.e. --repeat-only). This utility
is also able to calculate the distance of each protein (CDS) to the nearest
repeat_region. If there is no repeat_region in the particular record, the
distance is set to numpy.Infinity. Additionally, this utility can fetch
predicted conserved protein domains of each protein by accessing the NCBI
Conserved Protein database via POST request.

Author: Kyung Hoi (Joseph) Min
Written for the Arbor Biotechnologies coding challenge.

-------------------------------------------------------------------------------
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
'''

import os
import time
import requests
import numpy as np
import pandas as pd
from Bio import SeqIO

class CDS:
    '''
    Protein class.
    This is a dummy class that holds useful information
    about a particular protein.
    '''
    # Columns to include in output string.
    COLUMNS = ['locus_tag',
               'protein_id',
               'record_id',
               'location',
               'dist_to_repeat',
               'product',
               'translation']

    def __init__(self, record_id, location):
        '''
        Constructor.
        '''
        # Every CDS feature in the gbff file HAS an associated record and
        # a location.
        self.record_id = record_id
        self.location = location

        # But the may not have the other ones.
        self.locus_tag = None
        self.protein_id = None
        self.translation = None
        self.product = None
        self.closest_repeat = None
        self.dist_to_repeat = None

    def get_formatted_data(self):
        '''
        Returns a list of the values of each element in the self.columns
        list.
        '''
        return [str(getattr(self, col)) for col in self.COLUMNS]

    def __str__(self):
        '''
        Returns a string representation of this object in csv format.
        '''
        out = ','.join(self.get_formatted_data())

        return out

    def __repr__(self):
        '''
        Returns a dictionary representation of this object.
        '''
        return str(self.__dict__)

def find_conserved_domains(cdss):
    '''
    Sends an HTTP POST request to the NCBI Conserved Domain Search
    (https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi) and fetches domain
    information about the protein.
    This function specifically uses the amino acid sequence of the cds.
    (cds.translation).
    If no translation is available, no request is sent.
    '''
    def generate_queries(cdss):
        '''
        Returns a generator of query strings from a list of CDS objects.
        We generate a multiple queries instead of one in the case that we
        are attempting to send more than 4000 proteins (in which case
        the request will be blocked by the server).
        For safety, this function splits the queries into blocks of maximum
        sequences.
        '''
        print('Generating query strings')

        # This query string will be in FASTA format.
        query = ''
        query_count = 0
        for cds in cdss:
            # We query only those CDS's with the translated sequence.
            if cds.translation:
                # Each query CDS starts with a header line, which starts with
                # the '>' character. We use locus_tag here because there are
                # some CDS's without a protein_id, but all CDS's seem to have a
                # locus_tag.
                header = '>{}'.format(cds.locus_tag)
                sequence = cds.translation

                # Add the data to the query string.
                # Make sure to make a new query if there is more than 2000
                # sequences already in this query.
                if query_count > 2000:
                    yield query + '\n'
                    query_count = 0
                    query = ''
                query += '{}\n{}\n'.format(header, sequence)
                query_count += 1

        yield query

    def send_post_and_get_id(base_url, args):
        '''
        Sends the arguments to the server and retrieves the search id (cdsid).
        '''
        print('Sending POST request to {}'.format(base_url))
        print('{} proteins'.format(args['queries'].count('>')))
        r = requests.post(base_url, data=args)

        # Extract cdsid (search id) for this request.
        cdsid = None
        for line in r.text.split('\n'):
            if line.startswith('#cdsid'):
                cdsid = line.strip().split('\t')[1]

        if not cdsid:
            raise Exception('Could not extract cdsid.')

        print('Sent search query {}'.format(cdsid))

        return cdsid

    def wait_for_result(base_url, cdsid):
        '''
        Given a search id (cdsid), waits until the server is done computing
        the conserved domains. When it is, returns the result.
        '''
        print('Waiting for {}'.format(cdsid))

        args = {'cdsid': cdsid,
                'tdata': 'hits'}

        # Check for completion every 5 seconds.
        while True:
            time.sleep(5)
            r = requests.post(base_url, data=args)

            # Extract status.
            look_for = '#status'
            status = int(r.text[r.text.find(look_for) + len(look_for) + 1])

            # Check status.
            if status == 0:
                print('\nCompleted {}'.format(cdsid))

                # Return only the domain info text, not the search status text.
                return r.text.split('\n\n')[1]

            elif status == -1:
                raise Exception('Failed to retreive status')
            elif status == 1:
                raise Exception('Invalid cdsid')
            elif status == 2:
                raise Exception('Invalid POST arguments')
            elif status == 3:
                print('.', end='', flush=True)
            elif status == 4:
                raise Exception('Queue Manager Service error')
            elif status == 5:
                raise Exception('Data corrupted or no longer available')
            else:
                raise Exception('Unrecognized status {}'.format(status))

    def parse_domains(text):
        '''
        Helper function to parse the domain information returned (as a str),
        into a Pandas.DataFrame object.
        '''
        cols = None
        data = []
        for line in text.strip().split('\n'):
            split = line.split('\t')

            # If we don't have the columns yet, extract.
            if not cols:
                # Rename first column 'Query' to 'locus_tag'
                cols = split
                cols[0] = 'locus_tag'
            else:
                # All the query cells will have the locus tag.
                split[0] = split[0].split('>')[1]
                data.append(split)

        return pd.DataFrame(data, columns=cols)

    print('Finding conserved domains of each protein')

    # Fixed URL to the Batch CD-Search server.
    base_url = 'https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi'

    # POST arguments.
    args = {'db'    : 'cdd',
            'smode' : 'auto',
            'useid1': 'true',
            'evalue': 0.01,
            'maxhit': 500,
            'tdata' : 'hits',
            'dmode' : 'std',  # rep, std, full
            'cddefl': 'false'}

    # list of cdsids.
    cdsids = []

    # Send all queries first, before waiting on them.
    for query in generate_queries(cdss):
        args['queries'] = query
        cdsid = send_post_and_get_id(base_url, args)
        cdsids.append(cdsid)
        break

    # Wait for all queries to finish.
    df = None   # Pandas.DataFrame for all results.
    for cdsid in cdsids:
        result = wait_for_result(base_url, cdsid)
        domains = parse_domains(result)

        if df is None:
            df = domains
        else:
            # Otherwise, append the result to the dataframe.
            df.append(domains)

    print('Successfully retrieved conserved domains')
    print('-' * 50)

    # Return the final dataframe.
    return df


def find_closest_repeat(cds, repeat_regions):
    '''
    Given a CDS object and a list of all repeat_regions in the record,
    finds the closest repeat region, along with the distance in nucleotides to
    that region. If the repeat_regions list is empty, the distance is
    calculated to be infinity.
    '''
    # Location of the cds.
    # Note that this is a Bio.SeqFeature.FeatureLocation object.
    # loc.start = starting position, loc.end = end position, loc.strand = 1/-1
    loc = cds.location
    start = loc.start
    end = loc.end
    strand = loc.strand

    # Filter out repeat regions on the other strand.
    filtered = repeat_regions
    # filtered = list(filter(lambda x: strand == x.location.strand,
                           # repeat_regions))

    # If the filtered list is empty, then the distance is just infinity.
    closest = None
    dist = np.Infinity
    if len(filtered) > 0:
        for region in filtered:
            region_start = region.location.start
            region_end = region.location.end

            # If the repeat region occurs after the cds.
            if region_start > end:
                current_dist = region_start - end

            # If the repeat region occurs before the cds.
            elif region_end < start:
                current_dist = start - region_end

            # Otherwise, it means that they overlap (is this possible?).
            else:
                current_dist = 0

            # Replace closest and dist with current region if it is closer.
            if dist > current_dist:
                closest = region
                dist = current_dist

        # Once we have looped through all the repeat regions,
        # we have found the closest one.
        cds.closest_repeat = closest
        cds.dist_to_repeat = dist

    else:
        cds.dist_to_repeat = np.Infinity


def parse_gbff(gbff, filt=[], repeat_only=False):
    '''
    Given a gbff file (str), parses the gbff file using Bio.SeqIO.parse.
    If the filter argument is given, only returns information about the
    protein that satisfies the filter.
    (i.e. proteins that have the filter as a substring in either its protein_id
    or locus_tag)
    '''
    print('Parsing {}'.format(gbff))

    # Parse gbff file with Bio.SeqIO.parse.
    all_cdss = []
    for record in SeqIO.parse(gbff, format='genbank'):
        record_id = record.id
        record_seq = record.seq  # Full sequence of this locus.

        print('Parsing record {}'.format(record_id))

        # Iterate through all the features and use only those with
        # type = CDS. We keep track of all the repeat regions we have
        # found so far to calculate the distance.
        cdss = []
        repeat_regions = []     # list of Bio.SeqFeature
        for feature in record.features:
            if feature.type == 'CDS':
                # We know that every CDS has a location.
                location = feature.location

                # Make a new CDS object.
                cds = CDS(record_id, location)

                # However, other features might be missing.
                # This is the list of fields to look for in the qualifiers.
                to_find = ['locus_tag', 'protein_id', 'translation', 'product']

                # For each of the terms in to_find, populate the corresponding
                # field in the CDS object if it exists in the feature.
                # Also, if a filter is given, check here if this cds satisfies
                # the filter.
                include = False
                for term in to_find:
                    if term in feature.qualifiers:
                        # The values in the qualifiers always seem to be a list
                        # of 1 element.
                        found = feature.qualifiers[term]
                        assert len(found) == 1
                        found = found[0]
                        setattr(cds, term, found)

                        if len(filt) > 0 and all(f in found for f in filt):
                            include = True

                # If a filter is given, add the cds to the list only if
                # it satisfies the filter.
                if len(filt) > 0:
                    if include:
                        cdss.append(cds)
                else:
                    cdss.append(cds)

            # When we encounter a repeat_region, keep track of its location
            # by appending it to the list.
            elif feature.type == 'repeat_region':
                repeat_regions.append(feature)

        # Once we have processed all the features in a particular record,
        # calculate the distance of each feature to the nearest repeat region.
        # Note that this is PER RECORD.
        for cds in cdss:
            find_closest_repeat(cds, repeat_regions)

        # Once finished, add all the cdss to the all_cdss list.
        all_cdss.extend(cdss)

        # If we only want sequences with distance < infinity to a repeat
        # region, filter one more time.
        if repeat_only:
            all_cdss = list(filter(lambda x: x.dist_to_repeat < np.Infinity,
                                   all_cdss))

    print('Successfully parsed {}'.format(gbff))
    print('-' * 50)

    return all_cdss

def cdss_to_df(cdss):
    '''
    Given a list of CDS objects, converts the list to a Pandas.DataFrame,
    with CDS.COLUMNS as the columns.
    '''
    columns = CDS.COLUMNS
    data = []
    for cds in cdss:
        data.append(cds.get_formatted_data())

    # Convert to dataframe.
    return pd.DataFrame(data, columns=columns)

def filter_by_domain(df_cdss, df_domains, filt=[]):
    '''
    Given a list of domain filters filt, extract only proteins (CDS) that
    satisfy the filter in the 'superfamily' or 'Short name' fields.
    '''
    to_check = ['Superfamily', 'Short name']

    # First, filter the domains dataframe.
    # bool_list is a list of booleans that determines whether or not to
    # include a particular row of df_domains.
    bool_list = np.array([True] * len(df_domains))
    for f in filt:
        # The specified filter can be in any one of the to_check columns,
        # so we need another boolean list to do an OR operation.
        bool_col = np.array([False] * len(df_domains))

        for col in to_check:
            bool_col = bool_col | df_domains[col].str.contains(f)

        bool_list = bool_list & bool_col

    # A list of locus tags that passed the filter.
    locus_tags = df_domains.locus_tag[bool_list].unique()

    # Return a tuple of cds and domain dataframes with proteins that passed
    # the filter.
    df_cdss_filtered = df_cdss[df_cdss.locus_tag.isin(locus_tags)]
    df_domains_filtered = df_domains[df_domains.locus_tag.isin(locus_tags)]

    return df_cdss_filtered, df_domains_filtered


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('gbff', type=str,
                        help='.gbff file containing the protein(s) of interest')
    parser.add_argument('-p', type=str, default='out',
                        help='prefix for the output file(s)')
    parser.add_argument('-g', '--gbff-filter', type=str, default='',
                        help='a comma-separated list of proteins to filter '
                        + 'for; can be protein_id, locus_tag, translation, '
                        + 'product')
    parser.add_argument('-d', '--domain-filter', type=str, default='',
                        help='a comma-separated list of strings to filter '
                        + 'proteins by their conserved domain; applied only '
                        + 'when domains are retrieved')
    parser.add_argument('-r', '--repeat-only', action='store_true',
                        help='only consider proteins that have < infinity '
                        + 'distance from a repeat_region')
    parser.add_argument('-o', '--override', action='store_true',
                        help='use this option to always fetch conserved '
                             + 'protein domains--MAY TAKE A LONG TIME')
    args = parser.parse_args()

    # Extract command-line inputs.
    gbff = args.gbff
    prefix = args.p
    gbff_filt = args.gbff_filter
    domain_filt = args.domain_filter
    repeat_only = args.repeat_only
    override = args.override

    # Make sure that gbff is a valid file.
    if not os.path.isfile(gbff):
        raise FileNotFoundError(gbff)

    # Convert comma-separated filter string to a list.
    gbff_filt_list = gbff_filt.split(',')

    # Filter out blank strings from the list because those can interfere with
    # correct filtering.
    gbff_filt_list = list(filter(lambda x: len(x) > 0, gbff_filt_list))

    # Do the same for the domain filter.
    domain_filt_list = domain_filt.split(',')
    domain_filt_list = list(filter(lambda x: len(x) > 0, domain_filt_list))

    # Parse gbff file.
    parsed = parse_gbff(gbff, filt=gbff_filt_list, repeat_only=repeat_only)

    # Find conserved domains of each CDS.
    # IMPORTANT: This step may take some time if there are a lot of sequences
    # to query. Use with caution if no filter was applied!!
    df_domains = None
    if len(gbff_filt_list) == 0 and not repeat_only and not override:
        # Don't retrieve conserved domains if there was no filter.
        print('No filter applied...please use the --override option, provide '
              + 'a filter with the --filter option or only consider proteins '
              + 'close to a repeat_region with the --repeat-only option.')
    else:
        # Fetch conserved domains.
        df_domains = find_conserved_domains(parsed)

    # Convert CDS objects to one DataFrame.
    df_cdss = cdss_to_df(parsed)

    # If we are given a domain filter, filter once more.
    if len(domain_filt_list) > 0 and df_domains is not None:
        df_cdss, df_domains = filter_by_domain(df_cdss, df_domains,
                                               domain_filt_list)

    # Output results.
    proteins_file = '{}_proteins.csv'.format(prefix)
    print('Writing proteins to {}'.format(proteins_file))
    df_cdss.to_csv(proteins_file, index=False)

    if df_domains is not None:
        domains_file = '{}_domains.csv'.format(prefix)
        print('Writing domains to {}'.format(domains_file))
        df_domains.to_csv(domains_file, index=False)

    print('DONE!\n')
