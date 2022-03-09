import math
import itertools
from tempfile import TemporaryDirectory
import subprocess
import os
import sys
import shutil
from functools import reduce

from Bio.Blast import NCBIXML, Record
import numpy as np
import igraph


class BlastStatsCalculator:
    def __init__(self, search_lambda, search_kappa, search_space, match_score, mismatch_penalty,
                 gap_existence_penalty, gap_extension_penalty):
        self.search_lambda = search_lambda
        self.ln_kappa = math.log(search_kappa)
        self.search_space = search_space
        self.match_score = abs(match_score)
        self.mismatch_penalty = -1 * abs(mismatch_penalty)
        self.gap_existence_penalty = -1 * abs(gap_existence_penalty)
        self.gap_extension_penalty = -1 * abs(gap_extension_penalty)
        self.ln2 = math.log(2)

    def calculate_alignment_score(self, num_identities, num_mismatched, num_gaps, num_gap_openings):
        return num_identities * self.match_score + \
               num_mismatched * self.mismatch_penalty + \
               num_gaps * self.gap_extension_penalty + \
               num_gap_openings * self.gap_existence_penalty

    def calculate_bit_score(self, num_identities, num_mismatched, num_gaps, num_gap_openings):
        S = self.calculate_alignment_score(num_identities, num_mismatched, num_gaps, num_gap_openings)
        return (self.search_lambda * S - self.ln_kappa) / self.ln2

    def calculate_bit_score_from_alignment_score(self, alignment_score):
        return (self.search_lambda * alignment_score - self.ln_kappa) / self.ln2

    def calculate_e_value(self, num_identities, num_mismatched, num_gaps, num_gap_openings):
        bit_score = self.calculate_bit_score(num_identities, num_mismatched, num_gaps, num_gap_openings)

        try:
            denom = 2 ** bit_score
        except OverflowError as e:
            return 0

        return self.search_space / denom

    def calculate_e_value_from_bit_score(self, bit_score):
        try:
            denom = 2 ** bit_score
        except OverflowError as e:
            return 0

        return self.search_space / denom


def count_gap_openings(seq):
    gaps = 0
    i = 0
    while i < len(seq):
        if seq[i] == '-':
            gaps += 1
            i += 1
            while seq[i] == '-':
                i += 1
        else:
            i += 1
    return gaps


def split_hsp(hsp, pos, blast_stats_calc):
    '''
    :param pos: Position in HSP to split at (relative to start of HSP).
    :return: Two new HSPs
    '''
    def trim_ends(q_seq, s_seq, reverse=True):
        '''
        Trim ends of sequences when ends do not help alignment
        '''
        num_bases_to_trim = 0
        direction = -1 if reverse else 1

        for bases in list(zip(q_seq, s_seq))[::direction]:
            if '-' in bases or len(set(bases)) == 2:
                num_bases_to_trim += 1
            else:
                break

        if num_bases_to_trim == 0:
            return q_seq, s_seq, 0
        elif reverse:
            return q_seq[:-1 * num_bases_to_trim], s_seq[:-1 * num_bases_to_trim], num_bases_to_trim
        else:
            return q_seq[num_bases_to_trim:], s_seq[num_bases_to_trim:], num_bases_to_trim

    def get_alignment_stats(q_seq, s_seq):
        midline = ''
        num_gap_bases = 0
        num_identities = 0
        num_mismatches = 0

        for i in range(len(q_seq)):
            if '-' == q_seq[i] or '-' == s_seq[i]:
                midline += ' '
                num_gap_bases += 1
            elif q_seq[i] != s_seq[i]:
                midline += ' '
                num_mismatches += 1
            else:
                midline += '|'
                num_identities += 1
        return midline, num_identities, num_mismatches, num_gap_bases

    # Empty biopython HSP objects
    hsp_A = Record.HSP()
    hsp_B = Record.HSP()

    # Get sequences for new HSPs
    query_A, sbjct_A, num_bases_trimmed_A = trim_ends(hsp.query[:pos], hsp.sbjct[:pos])
    hsp_A.query = query_A
    hsp_A.sbjct = sbjct_A
    query_B, sbjct_B, num_bases_trimmed_B = trim_ends(hsp.query[pos:], hsp.sbjct[pos:], False)
    hsp_B.query = query_B
    hsp_B.sbjct = sbjct_B

    hsp_A.align_length = len(query_A)
    hsp_B.align_length = len(query_B)

    # Get alignment stats
    midline_A, num_identities_A, num_mismatches_A, num_gap_bases_A = get_alignment_stats(query_A, sbjct_A)
    midline_B, num_identities_B, num_mismatches_B, num_gap_bases_B = get_alignment_stats(query_B, sbjct_B)

    hsp_A.match = midline_A
    hsp_B.match = midline_B

    hsp_A.identities = num_identities_A
    hsp_A.positives = num_identities_A
    hsp_B.identities = num_identities_B
    hsp_B.positives = num_identities_B

    hsp_A.gaps = num_gap_bases_A
    hsp_B.gaps = num_gap_bases_B


    # query in positive frame?
    if hsp.frame[0] > 0:
        hsp_A.query_start = hsp.query_start
        hsp_A.query_end = hsp_A.query_start + (hsp_A.align_length - hsp_A.query.count('-') - 1)

        hsp_B.query_end = hsp.query_end
        hsp_B.query_start = hsp.query_end - (hsp_B.align_length - hsp_B.query.count('-') - 1)

    else:
        raise RuntimeError('Negative frame query?? We have not prepared for this situation!!!')

    # subject in positive frame?
    if hsp.frame[1] > 0:
        hsp_A.sbjct_start = hsp.sbjct_start
        hsp_B.sbjct_end = hsp.sbjct_end

        hsp_A.sbjct_end = hsp_A.sbjct_start + (hsp_A.align_length - hsp_A.sbjct.count('-') - 1)
        hsp_B.sbjct_start = hsp.sbjct_end - (hsp_B.align_length - hsp_B.sbjct.count('-') - 1)

    else:
        hsp_A.sbjct_start = hsp.sbjct_start
        hsp_A.sbjct_end = hsp_A.sbjct_start - (hsp_A.align_length - hsp_A.sbjct.count('-') - 1)

        hsp_B.sbjct_end = hsp.sbjct_end
        hsp_B.sbjct_start = hsp.sbjct_end + (hsp_B.align_length - hsp_B.sbjct.count('-') - 1)

    # Find number of gap openings
    openings_A = count_gap_openings(hsp_A.query) + count_gap_openings(hsp_A.sbjct)
    openings_B = count_gap_openings(hsp_B.query) + count_gap_openings(hsp_B.sbjct)

    # Now calculate score, bitscore, and e-value
    hsp_A.score = blast_stats_calc.calculate_alignment_score(hsp_A.identities, num_mismatches_A, hsp_A.gaps, openings_A)
    hsp_B.score = blast_stats_calc.calculate_alignment_score(hsp_B.identities, num_mismatches_B, hsp_B.gaps, openings_B)
    hsp_A.bits = blast_stats_calc.calculate_bit_score_from_alignment_score(hsp_A.score)
    hsp_B.bits = blast_stats_calc.calculate_bit_score_from_alignment_score(hsp_B.score)
    hsp_A.expect = blast_stats_calc.calculate_e_value_from_bit_score(hsp_A.bits)
    hsp_B.expect = blast_stats_calc.calculate_e_value_from_bit_score(hsp_B.bits)

    hsp_A.frame = hsp.frame
    hsp_B.frame = hsp.frame

    if hsp_A.align_length == 0:
        hsp_A = None

    if hsp_B.align_length == 0:
        hsp_B = None

    return hsp_A, hsp_B


def exact_to_rel_pos(hsp, genome_pos, query=True):
    if query:
        start, end = sorted([hsp.query_start, hsp.query_end])
        seq = hsp.query
    else:
        start, end = sorted([hsp.sbjct_start, hsp.sbjct_end])
        seq = hsp.sbjct

    # out of bounds
    if genome_pos < start or genome_pos >= end:
        return -1

    if query or hsp.frame[1] > 0:
        i = 0
        while start <= genome_pos:
            if seq[i] != '-':
                start += 1

            i += 1
    else:
        i = 0
        while end > genome_pos:
            if seq[i] != '-':
                end -= 1

            i += 1

    return i


class Blast_Alignment_Graph():
    '''
    A graph in which nodes are represent HSPs of a specific BLAST alignment and an edge between HSPs signifies that both
    HSPs refer to non-overlapping regions in the query and subject sequence.
    '''
    def __init__(self, blast_alignment, blast_stats_calc, split_hsps=True, limit_at=2000, query_identifier='QUERY'):
        '''
        :param blast_alignment: An 'Alignment' object from parsing blast results with BioPython
        '''
        self.query_id = query_identifier
        self.subject_id = blast_alignment.hit_def
        self.alignment = blast_alignment
        self.blast_stats_calc = blast_stats_calc

        if split_hsps:
            self._overlap_split_init(limit_at)
        else:
            self._quick_init()

    def _quick_init(self):
        '''
        Initialize quickly without splitting HSPs at overlaps.
        Less accurate but faster than overlap_split_init.
        '''
        num_hsps = len(self.alignment.hsps)
        # First, we build two arrays: all start and stop locations for query, all start and stop locations for subject
        query_loc_array = []
        subject_loc_array = []
        for i in range(num_hsps):
            # Tuple = (HSP_id, location)
            query_loc_array.append((i, self.alignment.hsps[i].query_start))
            query_loc_array.append((i, self.alignment.hsps[i].query_end))
            subject_loc_array.append((i, self.alignment.hsps[i].sbjct_start))
            subject_loc_array.append((i, self.alignment.hsps[i].sbjct_end))

        # sort by location
        query_loc_array.sort(key=lambda tup: tup[1])
        subject_loc_array.sort(key=lambda tup: tup[1])

        # Determine Overlaps  TODO: Resolve issue of ties
        in_hit_query = [False] * num_hsps
        in_hit_subject = [False] * num_hsps
        adjacency_mat = np.matrix(np.zeros((num_hsps, num_hsps), np.bool))

        def mark_overlaps(in_hit_array, hsp_id):
            # At the start of HSP region
            if not in_hit_array[hsp_id]:
                # If start of this HSP is within other HSPs, record overlaps
                overlapping_elem = np.where(in_hit_array)[0]
                # Mark overlaps in adjacency matrix
                for ol_id in overlapping_elem:
                    adjacency_mat[hsp_id, ol_id] = True
                    adjacency_mat[ol_id, hsp_id] = True
                # Open HSP
                in_hit_array[hsp_id] = True
            # At the end of HSP region
            else:
                # Close HSP
                in_hit_array[hsp_id] = False

        for i in range(num_hsps * 2):
            q_hsp_id = query_loc_array[i][0]
            s_hsp_id = subject_loc_array[i][0]

            mark_overlaps(in_hit_query, q_hsp_id)
            mark_overlaps(in_hit_subject, s_hsp_id)

        # Flip matrix elements so edges represent a non-overlapping relationship between HSPs
        adjacency_mat = (adjacency_mat == False)
        np.fill_diagonal(adjacency_mat, False)

        # Create graph
        self.G = self.igraph_from_numpy_matrix(adjacency_mat)

    def _overlap_split_init(self, limit_at):
        '''
        Initilaze with spliting at overlaps. Use this method of initialization when calculating percent identity
        between two sequences.
        '''
        self._remove_interior_hsps()
        self._split_hsps_at_overlaps(on_query=True)
        self._split_hsps_at_overlaps(on_query=False)
        self._make_graph(limit_at)

    def _make_graph(self, limit_at):
        # Keep only the top XXXX HSPs
        top = []
        for i, hsp in enumerate(sorted(self.alignment.hsps, reverse=True, key=lambda x: x.bits)):
            if i == limit_at:
                break
            top.append(hsp)
        self.alignment.hsps = top

        # Determine Overlaps
        num_hsps = len(self.alignment.hsps)
        adjacency_mat = np.matrix(np.zeros((num_hsps, num_hsps), np.bool))

        def hsps_overlapping(hsp_A, hsp_B):
            '''
            Return true if given hsps are overlapping, false otherwise
            '''
            A1, A2 = sorted([hsp_A.query_start, hsp_A.query_end])
            B1, B2 = sorted([hsp_B.query_start, hsp_B.query_end])

            if A1 <= B2 and B1 <= A2:
                return True
            else:
                A1, A2 = sorted([hsp_A.sbjct_start, hsp_A.sbjct_end])
                B1, B2 = sorted([hsp_B.sbjct_start, hsp_B.sbjct_end])

                if A1 <= B2 and B1 <= A2:
                    return True
                else:
                    return False

        def nodes_with_same_neighbors(hsp_id_1, hsp_id_2):
            return np.array_equal(adjacency_mat[hsp_id_1], adjacency_mat[hsp_id_2])

        for hsp_id_1, hsp_id_2 in itertools.combinations(range(num_hsps), 2):
            if hsps_overlapping(self.alignment.hsps[hsp_id_1], self.alignment.hsps[hsp_id_2]):
                adjacency_mat[hsp_id_1, hsp_id_2] = True
                adjacency_mat[hsp_id_2, hsp_id_1] = True

        # Flip matrix elements so edges represent a non-overlapping relationship between HSPs
        adjacency_mat = (adjacency_mat == False)
        np.fill_diagonal(adjacency_mat, False)

        # Create graph
        self.G = self.igraph_from_numpy_matrix(adjacency_mat)

    def _split_hsps_at_overlaps(self, on_query=True):
        '''
        Split hsps based on overlap for either query or subject.
        (Must be done for both to have fully split set of HSPs)
        :param on_query: If true, split on query, otherwise split on subject
        '''
        num_hsps = len(self.alignment.hsps)
        loc_dict = {}
        for i in range(num_hsps):
            if on_query:
                locs = [self.alignment.hsps[i].query_start - 1, self.alignment.hsps[i].query_end]
            else:
                start, end = sorted([self.alignment.hsps[i].sbjct_start, self.alignment.hsps[i].sbjct_end])
                locs = [start - 1, end]
            for loc in locs:
                if loc in loc_dict.keys():
                    loc_dict[loc].append(i)
                else:
                    loc_dict[loc] = [i]

        open_hsps = [False] * num_hsps
        hsp_id_overwrites = {}
        for loc in sorted(loc_dict):
            # Open closed HSPs at this location
            hsps_to_open = []
            for hsp_id in loc_dict[loc]:

                while hsp_id in hsp_id_overwrites.keys():
                    hsp_id = hsp_id_overwrites[hsp_id]

                if not open_hsps[hsp_id] and self.alignment.hsps[hsp_id].align_length > 1:
                    hsps_to_open.append(hsp_id)
                else:
                    open_hsps[hsp_id] = False

            # Split HSPs on overlap
            overlapping_hsp_ids = np.where(open_hsps)[0]
            for overlapping_id in overlapping_hsp_ids:
                while overlapping_id in hsp_id_overwrites.keys():
                    overlapping_id = hsp_id_overwrites[overlapping_id]

                # Find split position of hsp
                split_loc = exact_to_rel_pos(self.alignment.hsps[overlapping_id], loc, query=on_query)

                # Split overlapping HSP to get two new HSPs
                new_hsp_L, new_hsp_R = split_hsp(self.alignment.hsps[overlapping_id], split_loc,
                                                 self.blast_stats_calc)
                # Close old HSP we just split
                open_hsps[overlapping_id] = False

                # HSPs which are on negative strand of subject need to be handled differently
                is_reverse_sub_hsp = not on_query and self.alignment.hsps[overlapping_id].frame[1] < 0
                if is_reverse_sub_hsp:
                    new_hsp_L, new_hsp_R = new_hsp_R, new_hsp_L

                if new_hsp_R is not None:

                    new_hsp_R_id = len(self.alignment.hsps)
                    not_trimmed = (on_query and new_hsp_R.query_start - 1 == loc) or (
                                not on_query and (new_hsp_R.sbjct_start - 1 == loc or new_hsp_R.sbjct_end - 1 == loc))
                    # If bases were trimmed on the hsp to the right,
                    # we need to add the hsp open location to the loc_dict
                    if not not_trimmed:
                        open_hsps.append(False)
                        if on_query:
                            start = new_hsp_R.query_start - 1
                            if start in loc_dict.keys():
                                loc_dict[start].append(new_hsp_R_id)
                            else:
                                loc_dict[start] = [new_hsp_R_id]

                        else:
                            if self.alignment.hsps[overlapping_id].frame[1] > 0:
                                start = new_hsp_R.sbjct_start - 1
                                if start in loc_dict.keys():
                                    loc_dict[start].append(new_hsp_R_id)
                                else:
                                    loc_dict[start] = [new_hsp_R_id]
                            else:
                                end = new_hsp_R.sbjct_end - 1
                                if end in loc_dict.keys():
                                    loc_dict[end].append(new_hsp_R_id)
                                else:
                                    loc_dict[end] = [new_hsp_R_id]
                    else:
                        # If no bases were trimmed
                        open_hsps.append(True)

                    hsp_id_overwrites[overlapping_id] = new_hsp_R_id
                    self.alignment.hsps.append(new_hsp_R)

                # We need to make sure we don't re-open the overlapping hsp when we get to it's end and
                # we don't consider the overlapping hsp at the end
                else:
                    hsp_id_overwrites[overlapping_id] = None

                    for x in loc_dict.keys():
                        if overlapping_id in loc_dict[x]:
                            loc_dict[x] = list(filter(lambda a: a != overlapping_id, loc_dict[x]))

                if new_hsp_L is not None:
                    open_hsps.append(False)
                    self.alignment.hsps.append(new_hsp_L)

            for id in hsps_to_open:
                open_hsps[id] = True

        # Overwrite the alignment's hsp array with split hsps
        new_hsps = []
        for hsp_id in set(np.where(self.alignment.hsps)[0]) - hsp_id_overwrites.keys():
            # TODO: E-Value check
            new_hsps.append(self.alignment.hsps[hsp_id])
        self.alignment.hsps = new_hsps

    def igraph_from_numpy_matrix(self, adjacency_mat):
        '''
        Create and return an igraph Graph based on a numpy array given as input.
        Not meant to work with all adjacency matrices, only undirected. Values in matrix should be True/False.
        '''
        G = igraph.Graph()

        # Confirm input is 2x2 matrix
        if len(adjacency_mat.shape) != 2:
            raise ValueError('Incorrect adjacency matrix shape.')

        # Confirm matrix is square
        if adjacency_mat.shape[0] != adjacency_mat.shape[1]:
            raise ValueError('Incorrect adjacency matrix shape.')

        G.add_vertices(adjacency_mat.shape[0])

        # Add weights to vertices
        for idx, v in enumerate(G.vs):
            v['weight'] = self.alignment.hsps[idx].identities

        G.add_edges(zip(*np.where(adjacency_mat)))

        return G

    def _remove_interior_hsps(self):
        '''
        Remove HSPs which fall completely within an other HSP and have a lower score than that HSP.
        '''
        def lies_within(a_start, a_end, b_start, b_end):
            '''
            :return: True if 'a' lies within 'b'
            '''
            return b_start <= a_start and b_end >= a_end

        hsps_to_delete = set()
        num_hsps = len(self.alignment.hsps)
        for hsp_id_1, hsp_id_2 in itertools.combinations(range(num_hsps), 2):
            hsp1 = self.alignment.hsps[hsp_id_1]
            hsp2 = self.alignment.hsps[hsp_id_2]
            hsp1_qstart, hsp1_qend = sorted([hsp1.query_end, hsp1.query_start])
            hsp2_qstart, hsp2_qend = sorted([hsp2.query_end, hsp2.query_start])
            hsp1_sstart, hsp1_send = sorted([hsp1.sbjct_end, hsp1.sbjct_start])
            hsp2_sstart, hsp2_send = sorted([hsp2.sbjct_end, hsp2.sbjct_start])
            if lies_within(hsp1_qstart, hsp1_qend, hsp2_qstart, hsp2_qend) and \
                    lies_within(hsp1_sstart, hsp1_send, hsp2_sstart, hsp2_send) and \
                    hsp1.bits < hsp2.bits:
                hsps_to_delete.add(hsp_id_1)
            elif lies_within(hsp2_qstart, hsp2_qend, hsp1_qstart, hsp1_qend) and \
                    lies_within(hsp2_sstart, hsp2_send, hsp1_sstart, hsp1_send) and \
                    hsp2.bits < hsp1.bits:
                hsps_to_delete.add(hsp_id_2)

        hsp_ids_to_keep = set(range(num_hsps)) - hsps_to_delete
        hsps_to_keep = [self.alignment.hsps[i] for i in hsp_ids_to_keep]
        self.alignment.hsps = hsps_to_keep

    def _remove_hsps(self, hsp_ids_to_remove, adjacency_mat):
        hsp_ids_to_keep = [x for x in range(len(self.alignment.hsps)) if x not in hsp_ids_to_remove]

        # Remove from HSP list
        self.alignment.hsps = [x for i, x in enumerate(self.alignment.hsps) if i not in hsp_ids_to_remove]

        # Remove from adjacency matrix
        new_matrix = adjacency_mat[hsp_ids_to_keep][:, hsp_ids_to_keep]
        return new_matrix

    def write_dimacs(self, filename):
        with open(filename, 'w') as fh:
            fh.write('p edge {} {}\n'.format(len(self.G.vs), len(self.G.es)))
            for idx, v in enumerate(self.G.vs):
                fh.write('n {} {}\n'.format(idx + 1, round(v['weight'])))

            for e in self.G.es:
                fh.write('e {} {}\n'.format(e.tuple[0] + 1, e.tuple[1] + 1))

    def write_metis(self, filename):
        with open(filename, 'w') as graph_fh, open('{}.weights'.format(filename), 'w') as weights_fh:
            graph_fh.write('{} {}\n'.format(len(self.G.vs), len(self.G.es)))
            for idx, v in enumerate(self.G.vs):
                graph_fh.write(' '.join([str(x.index + 1) for x in set(v.neighbors())]) + '\n')
                weights_fh.write('{} {}\n'.format(v.index, v['weight']))

    def get_max_nobs_pls(self, return_clique=False):
        '''
        Calculate the maximum non-overlapping bitscore sum (nobs). This will find all maximal cliques in the blast graph
        and calculate the total bit score for each. The highest bit score found is returned.

        return_clique: Return the HSPs in the clique w/ the max NOBS instead of the score.
        '''
        qx_loc = '/mnt/genomics/common_projects/software/qualex-ms-1.2/qualex-ms'
        atb_loc = '/mnt/genomics/common_projects/software/qualex-ms-1.2/converter/asc2bin'

        tmp = TemporaryDirectory()
        # tmp.name = '/state/partition1/AMD-Phage-Pipeline/exempliphi/'
        # ascii_graph_out = os.path.join(tmp.name, '%s.ag' % self.alignment.accession)
        ascii_graph_out = os.path.join(tmp.name, '%s.graph' % self.alignment.accession)
        self.write_metis(ascii_graph_out)
        # self.write_dimacs(ascii_graph_out)
        # bin_graph_out = os.path.join(tmp.name, '%s.bg' % self.alignment.accession)
        # subprocess.run([atb_loc, ascii_graph_out, bin_graph_out], check=True, stdout=subprocess.PIPE)
        # subprocess.run([qx_loc, bin_graph_out], stdout=subprocess.PIPE)

        # clique_idxs = []
        # with open("%s.sol" % bin_graph_out) as fh:
        #     for line in fh:
        #         if line.startswith('v'):
        #             clique_idxs.append(int(line.split()[1]))

        pls_loc = '/mnt/genomics/common_projects/software/open-pls-1.0/bin/pls'
        result = subprocess.run([
            pls_loc,
            '--algorithm=mwis',
            '--input-file=%s' % ascii_graph_out,
            '--weighted',
            '--use-weight-file'
        ], stdout=subprocess.PIPE)

        clique_idxs = None
        for line in result.stdout.decode('ascii').split('\n'):
            ls = line.split()
            if ls[0] == 'best-solution':
                clique_idxs = [int(x) for x in ls[2:]]
                break

        if return_clique:
            return [self.alignment.hsps[i] for i in clique_idxs]
        else:
            return sum([self.alignment.hsps[i].bits for i in clique_idxs])

    def get_max_nobs_cl(self, return_clique=False):
        CLIQUER_LOC = '/mnt/genomics/common_projects/software/cliquer-1.21/cl'

        tmp = TemporaryDirectory()
        ascii_graph_out = os.path.join(tmp.name, '%s.ag' % self.alignment.accession)
        self.write_dimacs(ascii_graph_out)

        try:
            result = subprocess.run(
                [CLIQUER_LOC, '-q', ascii_graph_out],
                check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60
            )

        except subprocess.TimeoutExpired as e:
            new_num_hsps = len(self.alignment.hsps) - 10
            self._make_graph(new_num_hsps)
            print('Blast alignment graph for {} not resolving in time. '
                  'Rerunning with {} HSPs.'.format(self.alignment.accession, new_num_hsps))
            shutil.move(ascii_graph_out, os.path.join('/mnt/genomics/mlueder/', '%s.ag' % self.alignment.accession))
            return self.get_max_nobs_cl(return_clique=return_clique)

        clique_idxs = [int(i) - 1 for i in result.stdout.decode().split(':')[1].split()]
        if return_clique:
            return [self.alignment.hsps[i] for i in clique_idxs]
        else:
            return sum([self.alignment.hsps[i].bits for i in clique_idxs])

    def get_max_nobs_tsmmwc(self, return_clique=False):
        TSM_MWC_LOC = '/mnt/genomics/common_projects/software/tsm-release/tsm-mwc'

        tmp = TemporaryDirectory()
        ascii_graph_out = os.path.join(tmp.name, '%s.wclq' % self.alignment.accession)
        self.write_dimacs(ascii_graph_out)

        try:
            result = subprocess.run(
                [TSM_MWC_LOC, ascii_graph_out],
                stdout=subprocess.PIPE, timeout=60
            )

        except subprocess.TimeoutExpired as e:
            new_num_hsps = len(self.alignment.hsps) - 10
            self._make_graph(new_num_hsps)
            print('Blast alignment graph for {} not resolving in time. '
                  'Rerunning with {} HSPs.'.format(self.alignment.accession, new_num_hsps))
            return self.get_max_nobs_tsmmwc(return_clique=return_clique)

        # shutil.move(ascii_graph_out, os.path.join('/mnt/genomics/mlueder/graphs',
        #                                           '{}_{}.ag'.format(self.query_id, self.subject_id)))

        for line in reversed(result.stdout.decode().splitlines()):
            if line.startswith('M'):
                clique_idxs = [int(x) - 1 for x in line.split()[1:]]

        # find highest scoring node (wlmc doesn't consider single nodes as maximal cliques)
        highest_num_idents = 0
        node_id = None
        for idx, hsp in enumerate(self.alignment.hsps):
            if hsp.identities > highest_num_idents:
                highest_num_idents = hsp.identities
                node_id = idx

        clique = [self.alignment.hsps[i] for i in clique_idxs]
        clique_total_ident = sum(hsp.identities for hsp in clique)
        if highest_num_idents >= clique_total_ident:
            clique_idxs = [node_id]
            clique = [self.alignment.hsps[node_id]]

        if return_clique:
            return clique
        else:
            return sum([self.alignment.hsps[i].bits for i in clique_idxs])

    def get_max_nobs_qx(self, return_clique=False):
        tmp = TemporaryDirectory()
        ascii_graph_out = os.path.join(tmp.name, '%s.ag' % self.alignment.accession)
        self.write_dimacs(ascii_graph_out)

        bin_graph_out = os.path.join(tmp.name, '%s.bg' % self.alignment.accession)
        subprocess.run(['asc2bin', ascii_graph_out, bin_graph_out], check=True, stdout=subprocess.PIPE)
        subprocess.run(['qualex-ms', bin_graph_out], stdout=subprocess.PIPE)

        clique_idxs = []
        with open("%s.sol" % bin_graph_out) as fh:
            for line in fh:
                if line.startswith('v'):
                    clique_idxs.append(int(line.split()[1]))

        if return_clique:
            return [self.alignment.hsps[i] for i in clique_idxs]
        else:
            return sum([self.alignment.hsps[i].bits for i in clique_idxs])

    def get_sequence_identity(self, query_len, quick=False, return_num_identities=False):
        '''
        Calculate sequence identity of query to subject sequence.
        '''
        # if len(self.alignment.hsps) > 120:
        if quick:
            hsps = self.get_max_nobs_qx(return_clique=True)
        else:

            hsps = self.get_max_nobs_tsmmwc(return_clique=True)

        total_ident = sum(hsp.identities for hsp in hsps)

        if return_num_identities:
            return total_ident

        return (total_ident / max([self.alignment.length, query_len])) * 100


class PairwiseProteinAlignment():
    '''
    Represents a pairwise alignment. We convert from needle or blast results to this before creating a
    Proteome_Comparison_Alignment_Graph.
    '''
    def __init__(self, protA_id, protB_id, protA_len, protB_len, score, num_ident, num_pos):
        self.protA_id, self.protB_id, self.protA_len, self.protB_len, self.score, self.num_ident, self.num_pos = \
            protA_id, protB_id, protA_len, protB_len, score, num_ident, num_pos
        self._calculate_hssp_score()

    @classmethod
    def create_from_BioMSA(cls, bioMSA):
        # def count_aa(seq):
        #     count = 0
        #     for pos in seq:
        #         if pos != '-':
        #             count += 1
        #     return count

        return cls(
            bioMSA._records[0].id,
            bioMSA._records[1].id,
            # count_aa(bioMSA._records[0].seq),
            # count_aa(bioMSA._records[1].seq),
            len(bioMSA._records[0].seq),
            len(bioMSA._records[1].seq),
            bioMSA.annotations['score'],
            bioMSA.annotations['identity'],
            bioMSA.annotations['similarity']
        )

    def _calculate_hssp_score(self):
        L = max([self.protA_len, self.protB_len])
        PID = (self.num_ident * 100 / L)

        if L <= 11:
            self.hssp_score = PID - 100
        elif L <= 450:
            self.hssp_score = PID - (480 * math.pow(L, -0.32 * (1 + math.exp(-L/1000))))
        else:
            self.hssp_score = PID - 19.5


class Proteome_Comparison_Alignment_Graph():
    def __init__(self, alignments):
        '''
        :param alignments: python list of PairwiseProteinAlignment.
        '''
        self.alignments = alignments

        self.phageA_accession = self.alignments[0].protA_id[:self.alignments[0].protA_id.find('_')]
        self.phageB_accession = self.alignments[0].protB_id[:self.alignments[0].protB_id.find('_')]

        num_alignments = len(self.alignments)
        adjacency_mat = np.matrix(np.zeros((num_alignments, num_alignments), np.bool))

        def alignments_mutually_exclusive(alignment_A, alignment_B):
            return len({alignment_A.protA_id, alignment_A.protB_id, alignment_B.protA_id, alignment_B.protB_id}) < 4

        for alignment_A_id, alignment_B_id in itertools.combinations(range(num_alignments), 2):
            if alignments_mutually_exclusive(self.alignments[alignment_A_id], self.alignments[alignment_B_id]):
                adjacency_mat[alignment_A_id, alignment_B_id] = True
                adjacency_mat[alignment_B_id, alignment_A_id] = True

        adjacency_mat = (adjacency_mat == False)
        np.fill_diagonal(adjacency_mat, False)

        # Create graph
        self.G = self.igraph_from_numpy_matrix(adjacency_mat)

    def igraph_from_numpy_matrix(self, adjacency_mat):
        '''
        Create and return an igraph Graph based on a numpy array given as input.
        Not meant to work with all adjacency matrices, only undirected. Values in matrix should be True/False.
        '''
        G = igraph.Graph()

        # Confirm input is 2x2 matrix
        if len(adjacency_mat.shape) != 2:
            raise ValueError('Incorrect adjacency matrix shape.')

        # Confirm matrix is square
        if adjacency_mat.shape[0] != adjacency_mat.shape[1]:
            raise ValueError('Incorrect adjacency matrix shape.')

        G.add_vertices(adjacency_mat.shape[0])

        # Add weights to vertices
        for idx, v in enumerate(G.vs):
            v['weight'] = self.alignments[idx].score

        G.add_edges(zip(*np.where(adjacency_mat)))

        return G

    def get_max_nobs_cl(self):
        CLIQUER_LOC = '/mnt/genomics/common_projects/software/cliquer-1.21/cl'

        tmp = TemporaryDirectory()
        tmp.name = '/mnt/phage_01/pipeline_output/classiphi'
        ascii_graph_out = os.path.join(tmp.name, '{}-{}.graph'.format(self.phageA_accession, self.phageB_accession))
        self.write_dimacs(ascii_graph_out)
        result = subprocess.run(
            [CLIQUER_LOC, '-q', ascii_graph_out],
            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        clique_idxs = [int(i) - 1 for i in result.stdout.decode().split(':')[1].split()]
        return [self.alignments[i] for i in clique_idxs]

    def get_max_nobs_pls(self):
        '''
        Calculate the maximum non-overlapping bitscore sum (nobs). This will find all maximal cliques in the blast graph
        and calculate the total bit score for each. The highest bit score found is returned.

        return_clique: Return the HSPs in the clique w/ the max NOBS instead of the score.
        '''
        tmp = TemporaryDirectory()
        # tmp.name = '/mnt/phage_01/pipeline_output/classiphi'
        # ascii_graph_out = os.path.join(tmp.name, '%s.ag' % self.alignment.accession)
        ascii_graph_out = os.path.join(tmp.name, 'protein_alignment.graph')
        self.write_metis(ascii_graph_out)

        pls_loc = '/mnt/genomics/common_projects/software/open-pls-1.0/bin/pls'
        result = subprocess.run([
            pls_loc,
            '--algorithm=mwis',
            '--input-file=%s' % ascii_graph_out,
            '--weighted',
            '--use-weight-file'
        ], stdout=subprocess.PIPE)

        clique_idxs = None
        for line in result.stdout.decode('ascii').split('\n'):
            ls = line.split()
            if ls[0] == 'best-solution':
                clique_idxs = [int(x) for x in ls[2:]]
                break

        return [self.alignments[i] for i in clique_idxs]

    def write_dimacs(self, filename):
        with open(filename, 'w') as fh:
            fh.write('p edge {} {}\n'.format(len(self.G.vs), len(self.G.es)))
            for idx, v in enumerate(self.G.vs):
                fh.write('n {} {}\n'.format(idx + 1, round(v['weight'])))

            for e in self.G.es:
                fh.write('e {} {}\n'.format(e.tuple[0] + 1, e.tuple[1] + 1))

    def write_metis(self, filename):
        with open(filename, 'w') as graph_fh, open('{}.weights'.format(filename), 'w') as weights_fh:
            graph_fh.write('{} {}\n'.format(len(self.G.vs), len(self.G.es)))
            for idx, v in enumerate(self.G.vs):
                graph_fh.write(' '.join([str(x.index + 1) for x in set(v.neighbors())]) + '\n')
                weights_fh.write('{} {}\n'.format(v.index, v['weight']))

    def get_stats(self):
        clique = self.get_max_nobs_cl()
        # Here we define homologous as having an hssp score above 2
        num_homologous = 0
        for pa in clique:
            if pa.hssp_score > 2:
                num_homologous += 1

        return num_homologous

        # OLD METHOD
        # # Total proteome identity
        # total_num_ident = sum([x.num_ident for x in clique])
        # total_max_prot_len = sum([max([x.protA_len, x.protB_len]) for x in clique])
        # proteome_identity = float(total_num_ident) / float(total_max_prot_len)
        #
        # # Number of homologous genes
        # num_homo_at_30_ident, num_homo_at_30_pos = 0, 0
        # for pairwise_alignment in clique:
        #     max_len = max([pairwise_alignment.protA_len, pairwise_alignment.protB_len])
        #     identity = float(pairwise_alignment.num_ident) / float(max_len)
        #     positivity = float(pairwise_alignment.num_pos) / float(max_len)
        #     if identity >= 0.3:
        #         num_homo_at_30_ident += 1
        #
        #     if positivity >= 0.3:
        #         num_homo_at_30_pos += 1
        #
        # return proteome_identity, num_homo_at_30_ident, num_homo_at_30_pos


def main():
    # bc = BlastStatsCalculator(search_lambda=0.625, search_kappa=0.41, search_space=78060977995, match_score=2,
    #                           mismatch_penalty=3, gap_existence_penalty=5, gap_extension_penalty=2)
    # print(bc.calculate_e_value(50, 6, 10, 2))

    # grice analysis
    # blast_records = NCBIXML.parse(open('/share/apps/grice/merrell_meta/blast_results/CP10700.blastn.xml'))
    # blast_records = NCBIXML.parse(open('/share/apps/grice/merrell_meta/blast_results/CPUSU1.blastn.xml'))
    # blast_records = NCBIXML.parse(open('/share/apps/grice/merrell_meta/blast_results/SA2014N.blastn.xml')) Num HSPs before limiting = 34566
    # blast_records = NCBIXML.parse(open('/share/apps/grice/merrell_meta/blast_results/SurvivorA.blastn.xml')) # Num HSPs before limiting = 34296
    # blast_records = NCBIXML.parse(open('/share/apps/grice/merrell_meta/blast_results/SurvivorB.blastn.xml')) # Num HSPs before limiting = 34578
    blast_records = NCBIXML.parse(open(sys.argv[1]))

    # blast_records = NCBIXML.parse(open('/mnt/genomics/mlueder/blastn_results.xml'))
    # blast_records = NCBIXML.parse(open('/Users/mlueder/Desktop/junk/fake_blast_results.xml'))

    subject_ids = []
    SPACER = 1e10
    full_alignment = None
    bs = None
    total_query_len = 0
    total_subject_len = 0

    for rec_i, record in enumerate(blast_records):
        # print('lambda = {}'.format(record.ka_params[0]))
        # print('kappa = {}'.format(record.ka_params[1]))
        # print('ss = {}'.format(record.effective_search_space))

        total_query_len += record.query_length

        bs = BlastStatsCalculator(
            search_lambda=record.ka_params[0],
            search_kappa=record.ka_params[1],
            search_space=record.effective_search_space,
            match_score=record.sc_match,
            mismatch_penalty=record.sc_mismatch,
            gap_existence_penalty=record.gap_penalties[0],
            gap_extension_penalty=record.gap_penalties[1]
        )

        for alignment in record.alignments:

            if alignment.hit_id not in subject_ids:
                subject_i = len(subject_ids)
                subject_ids.append(alignment.hit_id)
                total_subject_len += alignment.length
            else:
                subject_i = subject_ids.index(alignment.hit_id)

            if subject_i == 0 and rec_i == 0:
                full_alignment = alignment

            else:
                for hsp in alignment.hsps:
                    hsp.sbjct_start = int(hsp.sbjct_start + SPACER * subject_i)
                    hsp.sbjct_end = int(hsp.sbjct_end + SPACER * subject_i)
                    hsp.query_start = int(hsp.query_start + SPACER * rec_i)
                    hsp.query_end = int(hsp.query_end + SPACER * rec_i)
                    full_alignment.hsps.append(hsp)

    full_alignment.length = total_subject_len
    g = Blast_Alignment_Graph(full_alignment, bs, limit_at=int(sys.argv[2]))
    print(g.get_sequence_identity(total_query_len))

    # for rec_i, record in enumerate(blast_records):
    #     bs = BlastStatsCalculator(
    #         search_lambda=record.ka_params[0],
    #         search_kappa=record.ka_params[1],
    #         search_space=record.effective_search_space,
    #         match_score=record.sc_match,
    #         mismatch_penalty=record.sc_mismatch,
    #         gap_existence_penalty=record.gap_penalties[0],
    #         gap_extension_penalty=record.gap_penalties[1]
    #     )
    #
    #     for alignment in record.alignments:
    #
    #         print('Query: %s' % record.query)
    #         print('Subject: %s' % alignment.title)
    #         g = Blast_Alignment_Graph(alignment, bs)
    #         print(g.get_sequence_identity(record.query_length))
    print('Complete...')


if __name__ == '__main__':
    main()