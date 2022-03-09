# Copyright 2018-2019 Leidos Inc. Naval Medical Research Center Biological Defense Research Directorate
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
from ete3 import Tree
from ete3 import TreeStyle, TextFace, NCBITaxa
import os
import sys
import unicodedata
import re
from Bio import SeqIO
import argparse
from Bio.Blast import NCBIXML
import sqlite3
import itertools
import shutil
from multiprocessing import Pool
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from exempliphi.pipeline import SeqRecord_Header_Handler
from blast_graph.blast_alignment_graph import Blast_Alignment_Graph


def getScore(alignment):
    '''
    Creates a blast alignment graph and returns max nobs.
    '''
    bg = Blast_Alignment_Graph(alignment, split_hsps=False, blast_stats_calc=None)
    # return bg.get_max_nobs_tsmmwc()
    return bg.get_max_nobs_qx()


class Bit_Weighted_Phylo_Tree():

    def __init__(self, ncbi):
        self.t = Tree()
        self.ncbi = ncbi
        self.num_hits = 0

    def add_new_lineage(self, lineage, bit_score):
        '''
        Add weighted lineage to tree
        :param lineage: lineage array as returned by NCBITaxa.get_lineage()
        :param bit_score: bit score of the hit
        '''
        def recursive_helper(parent, lineage, ranks):
            # child_rank = self.ncbi.get_rank([lineage[0]])[lineage[0]]
            child_rank = ranks[lineage[0]]

            if lineage[0] in [x.name for x in parent.get_children()]:
                child = parent.search_nodes(name=lineage[0])[0]

                if child_rank != 'species':
                    child.score = child.score + bit_score

                else:
                    # When we hit a species more than once, we update tree only if bit score is greater than previous hit
                    if child.score < bit_score:
                        self.detach(child.name)
                        child = parent.add_child(name=lineage[0])
                        child.add_feature('score', bit_score)

            else:
                child = parent.add_child(name=lineage[0])
                child.add_feature('score', bit_score)
                if len(lineage) == 1 or child_rank == 'species':
                    self.num_hits += 1

            # Do not go past species rank
            if len(lineage) > 1 and child_rank != 'species':
                recursive_helper(child, lineage[1:], ranks)

        recursive_helper(self.t, lineage, self.ncbi.get_rank(lineage))

    def classify(self, alpha):
        '''
        Descends through tree, following the highest scoring nodes, ignoring at most 'alpha'% of the total bitscore
        :param alpha: Total amount of bitscore that can be ignored to achieve a more detailed classification
        :return: taxid (int) of determined classification

        The problems with this method:
        - Classification will be heavily biased towards what is in Genbank. Many matches to
          an organism only slightly related can outweigh one perfect match
        - If only one low quality hit from one organism, then it will classify over-confidently
        '''
        # If root node is missing, then no lineages have been added - return taxid for 'unclassified sequences' (12908)
        if len(self.t.search_nodes(name=1)) == 0:
            return 12908

        # Get total score from root
        total_score = self.t.search_nodes(name=1)[0].score

        '''
        If a tree has both virus (taxid=10239) and bacteria (taxid=2) the organism could be a prophage. If it is a 
        prophage, we want to ignore the bacterial part of the tree. To determine if it is a prophage, we set a cut-off 
        percentage for the bitscore needed in the viral part of the tree. If this cut-off is met, the viral part of the 
        tree is trimmed.
        '''
        viral_percent_cutoff = .05  # range = 0 - 1
        if len(self.t.search_nodes(name=10239)) > 0 and len(self.t.search_nodes(name=2)) > 0:
            if self.t.search_nodes(name=10239)[0].score / total_score > viral_percent_cutoff:
                self.detach(2)
                # Recalculate total score
                total_score = self.t.search_nodes(name=1)[0].score

        def recursive_helper(node):
            if node.is_leaf():
                return int(node.name)

            else:
                next = max(node.get_children(), key=lambda node: node.score)
                if next.score / total_score < 1 - alpha:
                    return int(node.name)

                else:
                    return recursive_helper(next)

        return recursive_helper(self.t)

    def detach(self, name):
        '''
        Special detach method which updates the score of parent nodes. If score of parent node goes to 0, it is removed.
        '''
        orig_node = self.t.search_nodes(name=name)[0]
        def recursive_helper(node):
            node.score = node.score - orig_node.score
            parent = node.up
            if node.score <= 0:
                node.detach()

            if not parent.is_root():
                recursive_helper(parent)

        recursive_helper(orig_node.up)


class Phylo_Tree_Drawer():
    # Ranks we always want to show
    sig_ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'subfamily', 'genus', 'species']

    def __init__(self, ncbi):
        self.t = Tree()
        self.style = TreeStyle()
        self.style.show_leaf_name = False
        self.style.show_scale = False
        self.ncbi = ncbi
        self.is_collapsed = False

    def add_new_lineage(self, lineage, num_reads):
        '''
        Add lineage to tree
        :param lineage: lineage array as returned by NCBITaxa.get_lineage()
        :param num_reads: the number of reads you are classifying as part of the lineage (additive)
        '''
        if self.is_collapsed:
            raise RuntimeError("You may not add new lineages after drawing or creating directories")

        def recursive_helper(parent, lineage):
            if lineage[0] in [x.name for x in parent.get_children()]:
                child = parent.search_nodes(name=lineage[0])[0]
            else:
                child = parent.add_child(name=lineage[0])

            if len(lineage) > 1:
                recursive_helper(child, lineage[1:])
            elif 'num_reads' not in child.features:
                child.add_feature('num_reads', num_reads)
            else:
                child.num_reads = child.num_reads + num_reads

        recursive_helper(self.t, lineage)

    def draw(self, full_tree_path, simplified_tree_path, significance_ratio):
        '''
        Output full tree image at full_tree_path.
        Output simplified tree image at simplified_tree_path.
        Significance_ratio is the percent abundance needed to be considered significant. (Significant nodes are
        highlighted and the is a simplified tree image to show significant nodes.)
        Only call after you are done adding lineages.
        Requires X server
        '''
        self.ncbi.annotate_tree(self.t)
        if not self.is_collapsed:
            self._collapse_tree()
            self.is_collapsed = True

        # calculate total # of reads
        total_reads = 0
        for node in self.t.traverse():
            if 'num_reads' in node.features:
                total_reads += node.num_reads

        sig_threshold = round(total_reads * significance_ratio)
        self._add_text_faces(sig_threshold, total_reads)

        # Draw full tree image
        self.t.render(full_tree_path, w=35, units='in', tree_style=self.style)

        # Draw simplified tree image
        self._collapse_tree(sig_threshold)
        self.t.render(simplified_tree_path, w=35, units='in', tree_style=self.style)

    def create_directories(self, path, sequence_dict):
        '''
        This function will create a directory tree in the shape of the phylogenetic tree and deposit sequences there
        :param path: Root directory path
        '''
        self.ncbi.annotate_tree(self.t)
        if not self.is_collapsed:
            self._collapse_tree()
            self.is_collapsed = True

        def create_folder(path):
            head, tail = os.path.split(path)
            if head and not os.path.isdir(head):
                create_folder(head)

            if not os.path.isdir(path):
                os.mkdir(path)

        def slugify(value):
            """
            Convert spaces to underscores. Remove characters that aren't alphanumerics, underscores, or hyphens.
            Convert to lowercase. Also strip leading and trailing whitespace.
            """
            value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore').decode('ascii')
            value = re.sub(r'[^\w\s-]', '', value).strip().lower()
            return re.sub(r'[-\s]+', '_', value)

        def recursive_helper(node, *folders):
            name_slug = slugify(node.sci_name)
            full_path = os.path.join(path, *folders, name_slug)
            create_folder(full_path)

            # If sequences belong to this taxon, deposit them here
            if node.taxid in sequence_dict:
                SeqIO.write(
                    sequence_dict[node.taxid], os.path.join(full_path, '%s_sequences.fasta' % name_slug), 'fasta'
                )

            for child in node.get_children():
                recursive_helper(child, *folders, name_slug)

        for top_level_node in self.t.get_children():
            recursive_helper(top_level_node)

    def _add_text_faces(self, highlight_treshold, total_reads):
        '''
        Add labels to the image
        :param highlight_treshold: If node has over this number of reads, highlight the textface
        '''
        for node in self.t.traverse():
            label_text = '%s\nRank: %s' % (node.sci_name, node.rank)
            num_reads = 0
            if 'num_reads' in node.features:
                num_reads = node.num_reads
                label_text = '%s\nNumber of Reads: %i' % (label_text, num_reads)

            if num_reads >= highlight_treshold:
                label_text = '%s\nAbundance: %s' % (label_text, '{:.2%}'.format(num_reads/total_reads))

            face = TextFace(label_text)

            if num_reads >= highlight_treshold:
                face.background.color = "Moccasin"

            node.add_face(face, column=0)

    def _collapse_tree(self, min_reads=1):
        '''
        Remove nodes which do not have at least 'min_reads'(int) assigned and are not a significant rank
        '''
        def recursive_helper(node):
            children = node.get_children()

            num_reads = 0
            if 'num_reads' in node.features:
                num_reads = node.num_reads

            has_sig_ancestor = False
            for child in children:
                if recursive_helper(child):
                    has_sig_ancestor = True

            if num_reads >= min_reads:
                return True

            if not (node.rank in self.sig_ranks and has_sig_ancestor):
                node.delete(prevent_nondicotomic=False)

            return has_sig_ancestor

        for child in self.t.get_children():
            recursive_helper(child)


def parse_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument(
        '--taxdb_loc',
        dest='taxdb_loc',
        required=True,
        help='Location of "taxa.sqlite". Include database filename and extension.',
        type=str
    )

    arg_parser.add_argument(
        '--blast_results_contigs',
        dest='blast_results_contigs',
        required=True,
        help='Path to xml output from blasting assembly contigs',
        type=str
    )

    arg_parser.add_argument(
        '--blast_results_reads',
        dest='blast_results_reads',
        required=True,
        help='Path to xml output from blasting unmapped reads',
        type=str
    )

    arg_parser.add_argument(
        '--output',
        dest='output_dir',
        required=True,
        help='Directory to output results to',
        type=str
    )

    arg_parser.add_argument(
        '--full_tree_name',
        dest='full_tree',
        required=True,
        help='Name of the full tree image file (including extension)',
        type=str
    )

    arg_parser.add_argument(
        '--simplified_tree_name',
        dest='simp_tree',
        required=True,
        help='Name of the simplified tree image file (including extension)',
        type=str
    )

    arg_parser.add_argument(
        '--sig_threshold',
        dest='sig_threshold',
        required=True,
        help='[0-1] Abundance of reads a taxonomic classification needs to be considered significant. (0.50 = 50%)',
        type=float
    )

    arg_parser.add_argument(
        '--contigs',
        dest='contigs',
        required=True,
        help='Path to fasta file containing assembled contigs',
        type=str
    )

    arg_parser.add_argument(
        '--unmapped_reads',
        dest='unmapped_reads',
        required=True,
        help='Path to fasta file containing reads which did not map to the assembly',
        type=str
    )

    arg_parser.add_argument(
        '--num_procs',
        dest='num_procs',
        default=1,
        help='Number of processors to use',
        type=int
    )

    return arg_parser.parse_args()

if __name__ == '__main__':
    print('Running taxonomy module.')

    args = parse_args()

    ncbi = NCBITaxa(dbfile=args.taxdb_loc)

    db_con = sqlite3.connect(args.taxdb_loc)
    c = db_con.cursor()

    # Find the best classification for each contig
    blast_records_contigs = NCBIXML.parse(open(args.blast_results_contigs))
    blast_records_reads = NCBIXML.parse(open(args.blast_results_reads))
    classifications = {}  # { contig_name: taxid }

    p = Pool(args.num_procs)

    print('Parsing BLAST results and classifying query sequences...')
    for record in itertools.chain(blast_records_contigs, blast_records_reads):
        tree = Bit_Weighted_Phylo_Tree(ncbi)

        alignment_scores = p.map(getScore, record.alignments)
        # alignment_scores = []
        # for align in record.alignments:
        #     alignment_scores.append(getScore(align))

        for i in range(len(record.alignments)):
            alignment = record.alignments[i]
            # If 'contig' or 'scaffold' is in name sequence is most likely misclassified, so ignore
            if 'contig' in alignment.hit_def.lower() or 'scaffold' in alignment.hit_def.lower():
                continue

            if alignment.accession[-2] == '.':
                alignment.accession = alignment.accession[:-2]

            result = c.execute('SELECT taxid FROM acc2taxid WHERE accession=?', (alignment.accession,))

            try:
                taxid = result.fetchone()[0]
                tree.add_new_lineage(ncbi.get_lineage(taxid), alignment_scores[i])
            except:
                print('Warning: Missing accession->taxid entry for id %s' % alignment.accession, file=sys.stderr)

        num_hit_to_alpha_map = {
            2: 0.3,
            3: 0.5,
            4: 0.55
        }
        alpha = num_hit_to_alpha_map.get(tree.num_hits, 0.65)
        taxid = tree.classify(alpha)
        classifications[record.query] = taxid

    # Remove old tree if present
    old_dir_path = os.path.join(args.output_dir, 'sequences')
    if os.path.isdir(old_dir_path):
        shutil.rmtree(old_dir_path)

    print('Creating output...')
    sequence_dict = {}  # { taxid: [SeqRecords] }
    tree_drawer = Phylo_Tree_Drawer(ncbi)
    for sequence in itertools.chain(SeqIO.parse(args.contigs, 'fasta'), SeqIO.parse(args.unmapped_reads, 'fasta')):
        if sequence.description in classifications:
            taxid = classifications[sequence.description]

            # Get number of reads which mapped to sequence
            h_handler = SeqRecord_Header_Handler(sequence)
            features = h_handler.get_header_features()

            # Contigs will have the number of reads in header
            if features and 'num_reads' in features:
                num_reads = int(features['num_reads'])
            else:
                num_reads = 1

            tree_drawer.add_new_lineage(ncbi.get_lineage(taxid), num_reads)

            if taxid in sequence_dict:
                sequence_dict[taxid].append(sequence)
            else:
                sequence_dict[taxid] = [sequence]

    tree_drawer.create_directories(os.path.join(args.output_dir, 'sequences'), sequence_dict)

    tree_drawer.draw(
        os.path.join(args.output_dir, args.full_tree),
        os.path.join(args.output_dir, args.simp_tree),
        args.sig_threshold
    )