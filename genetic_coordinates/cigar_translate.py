"""
Description:
A tool for translating transcript coordinates to reference coordinates using the CIGAR string.

Contributors:
20210119 - Larry Clos II (drlclos@gmail.com)

Assumptions:
- For a given genetic sequence (transcript) the sequence coordinate and index are the same and start at zero (0).

- CIGAR defines the the order and size of transcript alignment events over the range of a padded reference string ,
  but certain events may not map to unpadded reference coordinates.
  - CIGAR operators are preceded by non-zero integers.
  - CIGAR strings never start with an Insertion or Deletion operator.

- There is a 1:1 coordinate/index mapping for reference and transcript sequence strings in (mis)match event regions
  (i.e. map between characters).

- There is not a 1:1 coordinate/index mapping for transcript regions not mapped to the unpadded reference (insertions)
  or for reference regions missing from the transcript (deletions).

- If queried, transcript insertions or reference deletions should map to the closest upstream (mis)match position
  (left-aligned) within the compared strand.

- A reference coordinate can be translated from a transcript coordinate query with the following formula:
   <ref start position> + <transcript coordinate query> + sum(<relative indel event sizes in prior coordinates>)

Usage:
 $python translate_coordinates.py <path to transcript alignments> <path to coordinate queries> <path to output file>

Test:
 $python -m doctest -v translate_coordinates.py

"""

import csv
import re
import sys

from collections import namedtuple
from typing import List, Tuple, Union

# namedtuple for input transcripts is lightweight, memory efficient alternative to creating full class
Transcript = namedtuple('Transcript', 'name chromosome ref_start cigar transcript_length index_map')


def cigar_to_operators(cigar: str) -> List[Tuple[int, str]]:
    """Isolate individual CIGAR operators and their numeric values.

    >>> cigar_to_operators('8M7D6M2I2M11D7M')
    [(8, 'M'), (7, 'D'), (6, 'M'), (2, 'I'), (2, 'M'), (11, 'D'), (7, 'M')]
    >>> cigar_to_operators('X8M7D6M2I2M11D7M12')
    [(8, 'M'), (7, 'D'), (6, 'M'), (2, 'I'), (2, 'M'), (11, 'D'), (7, 'M')]
    >>> cigar_to_operators('8M7D6M2I2M11D7M12Q')
    [(8, 'M'), (7, 'D'), (6, 'M'), (2, 'I'), (2, 'M'), (11, 'D'), (7, 'M')]

    Args:
        cigar: CIGAR string representing relative alignment of a transcript to a reference.

    Returns:
        A list of tuples, each being (<size of alignment type (integer)>, <CIGAR operator code (string)>).

    """
    re_cigar = list((int(c[:-1]), c[-1]) for c in re.findall(r'[0-9]+[MIDNSHP=X]', cigar))
    return re_cigar


def get_length(cigar: str) -> int:
    """Deduce the length of a transcript from its CIGAR string.

    >>> get_length('2M2X4M7D6M2I2M11D7M')
    25
    >>> get_length('2M2X4M7D6M2I2M11D7M12Q')
    25

    Args:
        cigar: CIGAR string representing relative alignment of a transcript to a reference.

    Returns:
        Integer value of the length of the transcript.

    """
    re_cigar = cigar_to_operators(cigar)
    length = sum(o[0] for o in re_cigar if o[1] not in ['D', 'N', 'H', 'P'])
    return length


def make_index_map(cigar: str) -> List[Tuple[int, int]]:
    """

    Transform CIGAR operators to an effective, efficient, and convenient data structure that maps the relative index of
    CIGAR alignment events that alter the 1:1 correspondence of transcript and reference coordinates (indexes).

    >>> make_index_map('8M7D6M2I2M11D7M')
    [(7, 7), (13, -2), (17, 11)]
    >>> make_index_map('2M2X4M7D6M2I2M11D7M')
    [(7, 7), (13, -2), (17, 11)]
    >>> make_index_map('2M2X2I2M')
    [(3, -2)]

    Design:
        - The index map assumes the transcript coordinate and index are the same and start at zero (0).
        - (Mis)matches don't alter index/coordinate map, hence need not be encoded, saving memory and computation.
        - Index map entries have the index of a left-aligned indel event 'a', and the associated event size 'b'
           such that:
           - any query > 'a' must add 'b' to the reference coordinate.
           - an insertion event has a negative 'b' to prevent advance of reference coordinate from transcript index.
           - sum of all 'b's for entries whose 'a' < query gives the total coordinate adjustment for indel events.
           - queries found in insertion regions need adjustment by amount of query index extending into region.

    Args:
        cigar: CIGAR string representing relative alignment of a transcript to a reference.

    Returns:
        A list of tuples, each being...
         (<index of position immediately prior to indel (left-aligned)>, <indel length a.k.a. coordinate adjust>).

    """
    re_cigar = cigar_to_operators(cigar)
    # TODO: can add model build for reference as well
    transcript_index_map = []
    ref_size = 0
    trans_size = 0
    for size, code in re_cigar:
        if code in ['M', '=', 'X']:
            ref_size += size
            trans_size += size
            continue
        if code in ['D', 'N']:
            ref_size += size
            transcript_index_map.append((trans_size - 1, size))
        elif code in ['I', 'S']:
            transcript_index_map.append((trans_size - 1, -size))
            trans_size += size
        # TODO: consider nuance from other operators

    return transcript_index_map


def translate_coordinate(coordinate_query: int,
                         index_map: List[Tuple[int, int]],
                         ref_start: int,
                         transcript_length: int) -> Union[int, None]:
    """Determine the reference coordinate (index) of input coordinate query from transcript index map and parameters.

    >>> [translate_coordinate(n, [(7, 7), (13, -2), (17, 11)], 3, 25) for n in [0, 13, 14, 15, 16, 24]]
    [3, 23, 23, 23, 24, 43]

    Nothing should be returned from below...
    >>> translate_coordinate(25, [(7, 7), (13, -2), (17, 11)], 3, 25)

    Args:
        coordinate_query:  Transcript position (index) to traslate to reference coordinate.
        index_map:         A list of tuples, each being...
                           (<index of position immediately prior to indel (left-aligned)>, <indel length>).
        ref_start:         Reference position that aligns with transcript index 0.
        transcript_length: Length of transcript.

    Returns:
        Reference coordinate, or None when query outside of transcript index range.

    """
    # Confirm query is possible with given transcript length
    if coordinate_query >= transcript_length or coordinate_query < 0:
        return None

    # Calculate reference coordinate for query
    mapped_adjusts = list(a for a in index_map if coordinate_query > a[0])
    ref_coord = sum([ref_start,
                     coordinate_query,
                     sum(a[1] for a in mapped_adjusts)])
    # Adjust insertions to return left-aligned reference coordinate
    if mapped_adjusts and mapped_adjusts[-1][1] < 0:
        last_adjust = mapped_adjusts[-1]
        diff = last_adjust[0] - coordinate_query - last_adjust[1]
        if diff > 0:
            ref_coord += diff
    return ref_coord


def main(transcripts: str, queries: str, output: str) -> None:
    """Translate input transcript coordinate queries to reference coordinates (all 0-based).

    Args:
        transcripts: Path to input TSV file containing the transcripts, with columns for transcript name, chromosome,
                     starting position, and CIGAR mapping string.
        queries:     Path to input TSV file with a set of queries, with columns for transcript name and transcript
                     coordinate to translate to reference coordinate.
        output:      Path to output TSV file to write query results, with columns for transcript name, transcript
                     coordinate, chromosome, and reference coordinate determined from query.

    """
    # Parse inputs and preprocessing (precompute reused parameters)
    transcript_alignments = {}
    with open(transcripts, 'r') as f:
        for t in list(csv.reader(f, delimiter='\t')):
            transcript_alignments[t[0]] = Transcript(name=t[0],
                                                     chromosome=t[1],
                                                     ref_start=int(t[2]),
                                                     cigar=t[3],
                                                     transcript_length=get_length(t[3]),
                                                     index_map=make_index_map(t[3]))
    with open(queries, 'r') as q:
        raw_queries = list(csv.reader(q, delimiter='\t'))

    # Process queries and
    with open(output, 'w', newline='') as w:
        query_output = csv.writer(w, delimiter='\t')
        for transcript_name, query in raw_queries:
            query = int(query)
            transcript = transcript_alignments[transcript_name]

            ref_coord = translate_coordinate(coordinate_query=query,
                                             index_map=transcript.index_map,
                                             ref_start=transcript.ref_start,
                                             transcript_length=transcript.transcript_length)
            if ref_coord is None:
                raise ValueError(f'Transcript {transcript_name} has no coordinate (index) {query}')

            query_output.writerow([transcript_name, query, transcript.chromosome, ref_coord])


if __name__ == '__main__':
    """Execute translate_coordinates module from CLI."""
    main(*sys.argv[1:])

    import doctest
    doctest.testmod()
