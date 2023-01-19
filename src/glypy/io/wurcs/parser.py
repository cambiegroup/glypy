import re

try:
    from urllib import unquote
except ImportError:
    from urllib.parse import unquote

from glypy.composition import Composition
from glypy.structure import glycan, link as _link, glycan_composition
from glypy.io.tree_builder_utils import try_int

from .node_type import NodeTypeSpec
from .utils import base52, WURCSFeatureNotSupported


class WURCSParser(object):
    def __init__(self, line, structure_class=glycan.Glycan):
        self.line = unquote(line)
        self.structure_class = structure_class
        self.version = None
        self.node_type_count = None
        self.node_count = None
        self.edge_count = None
        self.node_type_map = {}
        self.node_index_to_node = {}
        self.glyph_to_node_index = {}
        self.has_uncertain_linkages = False

    def parse_node_type_section(self, section):
        # Extract UniqueRES list
        # a2122h-1x_1-5_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5_2*NCC/3=O][a2112h-1b_1-5][a1221m-1a_1-5
        # Split by "]["
        node_types = section.split("][")
        for i, node_type in enumerate(node_types, 1):
            self.node_type_map[i] = NodeTypeSpec.parse(node_type, self.version)
        return self.node_type_map

    def parse_node_index_to_type_section(self, section):
        for i, index in enumerate(map(int, section.split('-'))):
            alpha = base52(i)
            mono = self.node_type_map[index].to_monosaccharide()
            mono.id = i
            self.node_index_to_node[i] = mono
            self.glyph_to_node_index[alpha] = i
        return self.node_index_to_node

    parse_connection = re.compile(r"([a-zA-Z]+)([0-9]+|\?)")

    def parse_connectivity_map(self, section):
        if not section:
            return False

        if "{" in section or "}" in section:
            links = section.split("_")
            if len(links) > 1 and len(set(links)) == 1:
                # This is a composition, everybody is ambiguously linked to everybody
                return False
            raise WURCSFeatureNotSupported("Braced Undefined Linkages are not supported")

        links = section.split("_")
        for link in links:
            has_ambiguity = "|" in link
            has_bridge = "*" in link
            if has_bridge:
                raise WURCSFeatureNotSupported("Bridging MAPs are not supported.")

            parent_link_def, child_link_def = link.split("-", 1)

            parent_spec = self.parse_connection.findall(parent_link_def)
            child_spec = self.parse_connection.findall(child_link_def)

            parent_glyph = parent_spec[0][0]
            child_glyph = child_spec[0][0]
            if self.glyph_to_node_index[child_glyph] < self.glyph_to_node_index[parent_glyph]:
                parent_spec, child_spec = child_spec, parent_spec

            if has_ambiguity:
                parent_nodes = []
                parent_positions = []
                for parent_glyph, parent_position in parent_spec:
                    parent_positions.append(try_int(parent_position) or -1)
                    parent_nodes.append(
                        self.node_index_to_node[self.glyph_to_node_index[parent_glyph]])

                child_nodes = []
                child_positions = []
                for child_glyph, child_position in child_spec:
                    child_positions.append(try_int(child_position) or -1)
                    child_nodes.append(
                        self.node_index_to_node[self.glyph_to_node_index[child_glyph]])
                bond = _link.AmbiguousLink(
                    parent_nodes, child_nodes, parent_position=parent_positions,
                    child_position=child_positions, parent_loss=Composition("H"),
                    child_loss=Composition("OH"))
                bond.find_open_position()
            else:
                parent_glyph, parent_position = parent_spec[0]
                child_glyph, child_position = child_spec[0]
                parent_position = try_int(parent_position) or -1
                child_position = try_int(child_position) or -1

                parent = self.node_index_to_node[self.glyph_to_node_index[parent_glyph]]
                child = self.node_index_to_node[self.glyph_to_node_index[child_glyph]]
                bond = _link.Link(
                    parent, child, parent_position=parent_position, child_position=child_position,
                    parent_loss=Composition("H"), child_loss=Composition("OH"))
        return True

    def _to_composition(self):
        gc = glycan_composition.GlycanComposition()
        for node in self.node_index_to_node.values():
            gc[glycan_composition.MonosaccharideResidue.from_monosaccharide(node)] += 1
        return gc

    def parse(self):
        # RegEx from WURCS Framework (java).
        # Permalink https://gitlab.com/glycoinfo/wurcsframework/-/blob/c2f50bd8c330d4d0d5b03c1a3d4d32169dd19ac2/src/main/java/org/glycoinfo/WURCSFramework/util/array/WURCSImporter.java#L62

        # String strExp = "WURCS=(.+)/(\\d+),(\\d+),(\\d+)(\\+)?/\\[(.+)\\]/([\\d\\-<>]+)/(.*)";

        # WURCS=2.0/7,10,9/[a2122h-1x_1-5_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5_2*NCC/3=O][a2112h-1b_1-5][a1221m-1a_1-5]/1-2-3-4-2-5-4-2-6-7/a4-b1_a6-j1_b4-c1_d2-e1_e4-f1_g2-h1_h4-i1_d1-c3|c6_g1-c3|c6
        # group(0)	WURCS=2.0/7,10,9/[a2122h-1x_1-5_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5_2*NCC/3=O][a2112h-1b_1-5][a1221m-1a_1-5]/1-2-3-4-2-5-4-2-6-7/a4-b1_a6-j1_b4-c1_d2-e1_e4-f1_g2-h1_h4-i1_d1-c3|c6_g1-c3|c6
        # group(1)	2.0
        # group(2)	7
        # group(3)	10
        # group(4)	9
        # group(5)	""
        # group(6)	a2122h-1x_1-5_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5_2*NCC/3=O][a2112h-1b_1-5][a1221m-1a_1-5
        # group(7)	1-2-3-4-2-5-4-2-6-7
        # group(8)	a4-b1_a6-j1_b4-c1_d2-e1_e4-f1_g2-h1_h4-i1_d1-c3|c6_g1-c3|c6

        pattern = re.compile(r"WURCS=(.+)/(\d+),(\d+),(\d+)(\+)?/\[(.+)]/([\d\-<>]+)/(.*)")
        wurcs_match = pattern.search(self.line)
        self.version = float(wurcs_match.group(1))
        self.node_type_count = wurcs_match.group(2)
        self.node_count = wurcs_match.group(3)
        self.edge_count = wurcs_match.group(4)
        self.has_uncertain_linkages = wurcs_match.group(5) == "+"
        self.parse_node_type_section(wurcs_match.group(6))
        self.parse_node_index_to_type_section(wurcs_match.group(7))

        if self.parse_connectivity_map(wurcs_match.group(8)):
            return self.structure_class(root=self.node_index_to_node[0], index_method='dfs', canonicalize=True)
        return self._to_composition()


def loads(text, structure_class=glycan.Glycan):
    """Parse a WURCS-encoded glycan structure from `text` into a :class:`~.Glycan`
    or :class:`~.GlycanComposition`.

    Parameters
    ----------
    text : str
        The WURCS string to parse
    structure_class : :class:`type`, optional
        The class to use to wrap the :class:`~.Monosaccharide` graph (the default is :class:`~.Glycan`)

    Returns
    -------
    :class:`~.Glycan` or :class:`~.GlycanComposition`
        The parsed result
    """
    parser = WURCSParser(text, structure_class=structure_class)
    structure = parser.parse()
    return structure
