"""
Description:
A tool for extracting and transforming HGVS annotations from REST API.

Contributors:
20210204 - Larry Clos II (drlclos@gmail.com)

Assumptions for prototype:
* Ensembl server used for HGVS REST API queries for human species.
* Faulty input HGVS should be identified and represented in the table.
* Duplicate input HGVS should be identified, but consolidated to a single line in output.
* Type hints are helpful for code interpretation.

Issues:
- Code is custom built for Ensembl REST API endpoint https://rest.ensembl.org/vep/human/hgvs, with little flexibility...
 - solution:
  - isolate API queries in Ensembl-specific functions.
 - followup questions:
  - what is value/cost in generalizing REST API access functions for many servers with unique REST fields and
    response data structures?
- employing POST method for HGVS list is convenient and minimizes API calls, but has no response for bad HGVS inputs...
 - solution:
  - employed GET method to obtain explicit error info for HGVS annotations not returned by POST.
 - followup questions:
  - is there value in determining different error modes and/or codes?
- constructing POST query params should be optimized in later dev.
- rest_server parameter passed around too much, with common URL params...
 -  solution:
  - make URL object at high level with common params, then pass it around instead.
- Ensembl annotation transform function should be made more robust.
 - 'keys_of_interest' could updated to a dictionary to map Ensembl response keys to internal standard names.
- no general configs

Install:
- execute the following to create the virtual environment:
 - $ cd <dir of this package cloned from github>
 - $ <path to python3.8> -m venv venv/3.8
 - for Linux:
  - $ source venv/3.8/bin/activate
 - for Windows:
  - $ source venv/3.8/Scripts/activate
 - $ pip install -U setuptools pip

"""

import argparse
import csv
from typing import Any, Dict, Iterator, List, Optional, Tuple

import requests


# Step 1 - Extract HGVS annotations from REST API
def extract_hgvs_annotations(input_hgvs: str, rest_server: str, headers: dict) -> dict:
    """Extract HGVS annotations from REST API.

    Notes:
        Combines duplicate HGVS input strings and sorts via set operation.
        Could be expanded with server-specific API function call logic and configs.

    Args:
        input_hgvs: Path to input text file with HGVS strings to query server.
        rest_server: URL of REST API server to extract HGVS annotations.
        headers: REST API headers

    Returns:
        Dictionary of HGVS-query keys with rest_server response JSON dictionary values.

    """
    with open(input_hgvs, 'r', newline='\n') as i:
        hgvs_strings = list(set(i.read().splitlines()))

    annotations = {hgvs: variant_annotation for hgvs, variant_annotation
                   in post_ensembl_hgvs_annotations(hgvs_input=hgvs_strings,
                                                    rest_server=rest_server,
                                                    headers=headers)}

    failed_post_queries = [hgvs for hgvs in hgvs_strings if hgvs not in annotations]
    for failed_hgvs in failed_post_queries:
        annotations[failed_hgvs] = get_ensembl_hgvs_annotations(hgvs_input=failed_hgvs,
                                                                rest_server=rest_server,
                                                                headers=headers)
    return annotations


def get_ensembl_hgvs_annotations(hgvs_input: str, rest_server: str, headers: dict) -> dict:
    """GET query Ensemble REST API on vep/human/hgvs endpoint with a single hgvs string.

    Notes:
        This isolates Ensembl-specific API GET interface logic.

    Args:
        hgvs_input: HGVS string to query server.
        rest_server: URL of REST API server to extract HGVS annotations.
        headers: REST API headers

    Returns:
        JSON dictionary returned from REST API GET query.

    """
    rest_server = rest_server + hgvs_input
    get_response = requests.get(rest_server, headers=headers)
    return get_response.json()


def post_ensembl_hgvs_annotations(hgvs_input: List[str], rest_server: str, headers: dict) -> Iterator[Tuple[str, dict]]:
    """POST query Ensemble REST API vep/human/hgvs endpoint with input hgvs strings.

    Notes:
        POST does not return annotations for faulty HGVS query (or other reason).
        GET used to obtain details of server error.
        This isolates Ensembl-specific API interface logic.

    Args:
        hgvs_input: List of HGVS strings to query server.
        rest_server: URL of REST API server to extract HGVS annotations.
        headers: REST API headers

    Yields:
        Valid HGVS query string and Ensembl annotation dictionary.

    """
    data = '{{"hgvs_notations": ["{0}"]}}'.format('", "'.join(hgvs_input))
    post_response = requests.post(rest_server, headers=headers, data=data)
    if post_response.status_code != 200:
        raise ValueError(post_response.text)
    for annotation in post_response.json():
        yield annotation['input'], annotation


# Step 2 - Transform REST API responses to internal data structure
def transform_ensembl_annotations(annotation: dict, keys_of_interest: List[str]) -> Dict[str, Any]:
    """Transform annotations returned from Ensembl REST API POST and GET calls.

    Notes:
          Can accommodate error responses -> returns 'NA' values for keys of interest.

    Args:
        annotation: JSON dictionary response from Ensembl REST API query.
        keys_of_interest:

    Returns:
        Dictionary of HGVS-query keys with transformed annotation dictionary values.

    """
    decoded_annotation = {}
    for key in keys_of_interest:
        # defer to top tier key non-container values, otherwise search for any nested key non-container values
        value = list(find_my_keys(key, annotation))
        if key in annotation and not isinstance(annotation[key], (list, dict)):
            decoded_annotation[key] = annotation[key]
        elif not value:
            decoded_annotation[key] = 'NA'
        else:
            decoded_annotation[key] = value
    return decoded_annotation


def find_my_keys(key: str, data: Any) -> Iterator[Any]:
    """Recursively find all values for instances of key in nested dictionary.

    Args:
        key: Key to search for in nested structure.
        data: Data structure to search for key.

    Yields:
        Non-dictionary values of any key:value pair found in data.

    """
    if isinstance(data, dict):
        for k, v in data.items():
            if not isinstance(v, (list, dict)) and k == key:
                yield v
            else:
                yield from find_my_keys(key, v)
    elif isinstance(data, list):
        for subelement in data:
            yield from find_my_keys(key, subelement)


# Step 3 - Generate outputs
def make_transformed_annotations_tsv(transformed_annotations: Dict[str, Dict[str, Any]],
                                     col_headers: list, output_tsv: str) -> None:
    """Generate flat output TSV file from transformed annotations.

    Notes:
        Keys of input annotations dict used for first column values ('query' header).
        This function could be more robust if server specific response keys were mapped in the prior transform step.
        Also could be more robust if separate function to make a transformed data frame (e.g. pandas).

    Args:
        transformed_annotations: Query HGVS string-keyed dictionary with transformed annotation dictionary values.
        col_headers: TSV column headers.
        output_tsv: Path to output TSV file

    """
    with open(output_tsv, 'w', newline='') as o:
        output_table = csv.writer(o, delimiter='\t')
        output_table.writerow(['query'] + col_headers)
        sorted_queries = sorted(hgvs for hgvs in transformed_annotations)
        for hgvs in sorted_queries:
            extract = transformed_annotations[hgvs]
            writable_output = [hgvs]
            for key in col_headers:
                value = extract[key]
                if isinstance(value, (str, int, float)):
                    writable_output.append(str(value))
                elif isinstance(value, list):
                    writable_output.append(','.join(sorted({str(found) for found in value})))
            output_table.writerow(writable_output)


# Module Workflow
def main(input_hgvs: str, rest_server: str, output_tsv: str) -> None:
    """Generate a flat TSV file from REST API HGVS-query annotation results.

    Args:
        input_hgvs: Path to input text file with HGVS strings to query server.
        rest_server: URL of REST API server to extract HGVS annotations.
        output_tsv: Path to output TSV file for transformed annotations.

    """
    # TODO: get keys of interest and headers from input config
    keys_of_interest = ['id',
                        'gene_id',
                        'gene_symbol',
                        'codons',
                        'impact',
                        'most_severe_consequence',
                        'clin_sig_allele',
                        'error']
    headers = {'content-type': 'application/json'}

    hgvs_annotations = extract_hgvs_annotations(input_hgvs,
                                                rest_server=rest_server,
                                                headers=headers)
    transformed = {hgvs: transform_ensembl_annotations(annotation=annotation,
                                                       keys_of_interest=keys_of_interest)
                   for hgvs, annotation in hgvs_annotations.items()}

    make_transformed_annotations_tsv(transformed_annotations=transformed,
                                     col_headers=keys_of_interest,
                                     output_tsv=output_tsv)


def cli(argv: Optional[List[str]] = None) -> None:
    """Command-line interface argument parser and module executor."""
    parser = argparse.ArgumentParser(prog='hgvs_annotations', description='ETL for HGVS annotations from Ensembl')
    parser.add_argument('-i', '--input',  help='Path to input text file with HGVS strings to query server.',
                        required=True)
    parser.add_argument('-o', '--output', help='Path to output TSV file for transformed annotations.',
                        required=True)
    parser.add_argument('-s', '--server', help='URL of REST API server to extract HGVS annotations.',
                        default='https://rest.ensembl.org/vep/human/hgvs')
    # TODO: add configuration file input to expand utility and customization of this module.

    args = parser.parse_args(argv)
    main(args.input, args.server, args.output)


if __name__ == '__main__':
    cli()
