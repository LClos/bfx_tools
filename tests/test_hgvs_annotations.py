# type: ignore
"""Tests for hgvs_annotations.py"""
from hgvs_annotations.hgvs_annotations import (find_my_keys,
                                               transform_ensembl_annotations,
                                               make_transformed_annotations_tsv,
                                               main)

from pathlib import Path
from pytest import mark

TEST_DATA = Path('tests')

TEST_DICT = {'A': 1,
             'B': [{'A': 2}, {'A': 3, 'B': 4}, {'C': 5}],
             'C': {'C': 6, 'D': '7'}
             }

TEST_ANNOTATIONS = {'A': {'a': 1, 'b': [2, 3, '4'], 'c': '5'},
                    'B': {'a': [6, '7'], 'b': 8, 'c': 9, 'd': '10'}
                    }


@mark.parametrize('key, data, expected',
                  [
                      ('A', TEST_DICT, [1, 2, 3]),
                      ('B', TEST_DICT, [4]),
                      ('C', TEST_DICT, [5, 6]),
                      ('D', TEST_DICT, ['7']),
                      ('X', TEST_DICT, [])
                  ]
                  )
def test_find_my_keys(key, data, expected):
    """Test find_my_keys function."""
    test_result = list(find_my_keys(key, data))
    assert test_result == expected


@mark.parametrize('data, keys, expected',
                  [
                      (TEST_DICT, ['A', 'B'], {'A': 1, 'B': [4]}),
                      (TEST_DICT, ['B', 'C'], {'B': [4], 'C': [5, 6]}),
                      (TEST_DICT, ['A', 'D', 'X'], {'A': 1, 'D': ['7'], 'X': 'NA'}),
                      ({}, ['A', 'D', 'X'], {'A': 'NA', 'D': 'NA', 'X': 'NA'})
                  ]
                  )
def test_transform_ensembl_annotations(data, keys, expected):
    """Test transform_ensembl_annotations function."""
    test_result = transform_ensembl_annotations(data, keys)
    assert test_result == expected


@mark.parametrize('data, col_headers, output_tsv, expected_file_result',
                  [
                      (TEST_ANNOTATIONS,
                       ['a', 'b', 'c'],
                       TEST_DATA / 'output/expected_simple_test_output.tsv',
                       TEST_DATA / 'input/expected_simple_test_output.tsv')
                  ]
                  )
def test_make_transformed_annotations_tsv(data, col_headers, output_tsv, expected_file_result):
    """Test make_transformed_annotations_tsv function."""
    make_transformed_annotations_tsv(data, col_headers, output_tsv)
    with open(output_tsv, 'r', newline='\n') as i:
        test_output_lines = [j.split('\t') for j in i.read().splitlines()]
    with open(expected_file_result, 'r', newline='\n') as i:
        expected_output_lines = [j.split('\t') for j in i.read().splitlines()]
    assert test_output_lines == expected_output_lines


@mark.parametrize('input_file, hgvs_server, output_file, expected_file_result',
                  [
                      (TEST_DATA / 'input/variants.txt',
                       'https://rest.ensembl.org/vep/human/hgvs',
                       TEST_DATA / 'output/full_test_result.tsv',
                       TEST_DATA / 'input/expected_full_test_output.tsv'
                       )
                  ]
                  )
def test_main(input_file, hgvs_server, output_file, expected_file_result):
    """Test make_transformed_annotations_tsv function."""
    main(input_file, hgvs_server, output_file)
    with open(output_file, 'r', newline='\n') as i:
        test_output_lines = [j.split('\t') for j in i.read().splitlines()]
    with open(expected_file_result, 'r', newline='\n') as i:
        expected_output_lines = [j.split('\t') for j in i.read().splitlines()]
    assert test_output_lines == expected_output_lines
