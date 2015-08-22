'''
A GlycoCT (Condensed) parser.
Supports RES, LIN, and un-nested REP sections.

.. code-block:: python

    >>>from glypy.io import glycoct
    >>>glycoct.loads("""RES
    1b:o-dman-HEX-0:0|1:aldi
    2b:a-lido-HEX-1:5|6:a
    3s:sulfate
    LIN
    1:1o(3+1)2d
    2:2o(2+1)3n
    3:2o(4+1)4d
    """).next()
    >>>
'''

import re
import logging
import warnings
from uuid import uuid4
from glypy.utils import opener, StringIO, enum, root as rootp
from glypy.utils.multimap import OrderedMultiMap
from glypy.structure import monosaccharide, substituent, link, glycan
from .format_constants_map import (anomer_map, superclass_map,
                                   link_replacement_composition_map, modification_map)
from glypy.composition import Composition

try:
    range = xrange
except:
    pass

logger = logging.getLogger("glycoct")

debug = False
import pdb


Glycan = glycan.Glycan
Monosaccharide = monosaccharide.Monosaccharide
Substituent = substituent.Substituent
Link = link.Link
AmbiguousLink = link.AmbiguousLink

START = "!START"
RES = "RES"
LIN = "LIN"
REP = "REP"
ALT = "ALT"
UND = "UND"

subsituent_start = "s"
base_start = "b"
repeat_start = "r"
alternative_start = "a"

#: Pattern for parsing the lines of the RES section corresponding
#: to individual |Monosaccharide| residues
res_pattern = re.compile(
    '''
    (?P<anomer>[abxo])?
    (?P<conf_stem>(?:-[dlx][a-z]+)+)?-?
    (?P<superclass>[A-Z]+)-?
    (?P<indices>[0-9x]+:[0-9x]+)
    (?P<modifications>(\|[0-9x]+:[0-9a-z]+)+)?
    ''', re.VERBOSE)

#: Pattern for parsing the potentially repeated |Configuration| and |Stem|
#: regions of the lines of the RES section.
conf_stem_pattern = re.compile(r'(?P<config>[dlx])(?P<stem>[a-z]+)')

#: Pattern for parsing modifications found on monosaccharide residue
#: lines in the RES section
modification_pattern = re.compile(r"\|?(\d+):([^\|;\n]+)")


#: Pattern for parsing |Link| lines found in the LIN section
link_pattern = re.compile(
    r'''(?P<doc_index>\d+)?:
    (?P<parent_residue_index>\d+)
    (?P<parent_atom_replaced>[odhnx])
    \((?P<parent_attachment_position>-?[0-9\-\|]+)[\+\-]
        (?P<child_attachment_position>-?[0-9\-\|]+)\)
    (?P<child_residue_index>\d+)
    (?P<child_atom_replaced>[odhnx])
        ''', re.VERBOSE)


#: Special truncation of the :data:`link_pattern` which is used on
#: REP header sections
internal_link_pattern = re.compile(
    r'''(?P<parent_residue_index>\d+)
    (?P<parent_atom_replaced>[odhnx])
    \((?P<parent_attachment_position>-?[0-9\-\|]+)[\+\-]
        (?P<child_attachment_position>-?[0-9\-\|]+)\)
    (?P<child_residue_index>\d+)
    (?P<child_atom_replaced>[odhnx])
    ''',
    re.VERBOSE)

#: Pattern for interpreting the REP# instance header section
rep_header_pattern = re.compile(
    r'''REP(?P<repeat_index>\d+):
    (?P<internal_linkage>.+)
    =(?P<lower_multitude>-?\d+)-(?P<higher_multitude>-?\d+)''', re.VERBOSE)

repeat_line_pattern = re.compile("^(?P<graph_index>\d+)r:r(?P<repeat_index>\d+)")


def parse_link(line):
    link_dict = link_pattern.search(line)
    if link_dict is not None:
        link_dict = link_dict.groupdict()
    else:
        raise GlycoCTError("Could not interpret link", line)
    id = link_dict['doc_index']
    parent_residue_index = link_dict['parent_residue_index']
    child_residue_index = link_dict['child_residue_index']

    parent_atom_replaced = link_replacement_composition_map[link_dict["parent_atom_replaced"]]
    parent_attachment_position = map(int, link_dict["parent_attachment_position"].split("|"))

    child_atom_replaced = link_replacement_composition_map[link_dict["child_atom_replaced"]]
    child_attachment_position = map(int, link_dict["child_attachment_position"].split("|"))

    return (id, parent_residue_index, parent_atom_replaced, parent_attachment_position,
            child_residue_index, child_atom_replaced, child_attachment_position)


def form_link(parent, child, parent_position, child_position, parent_loss, child_loss, id=None):
    #logger.debug("form_link %s", (parent_loss, child_loss, id))
    if parent.node_type is Substituent.node_type and\
     child.node_type is Monosaccharide.node_type:
        warnings.warn(
            "A monosaccharide with a substituent parent has been detected. "
            "These structures are not fully supported and may not traverse as expected "
            "by default.", stacklevel=7)
    if len(parent_position) > 1 or len(child_position) > 1:
        ambilink = AmbiguousLink(parent, child,
                                 parent_position=parent_position, child_position=child_position,
                                 parent_loss=parent_loss, child_loss=child_loss, id=id)
        ambilink.find_open_position()
    else:
        Link(parent, child, parent_position=parent_position[0],
             child_position=child_position[0],
             parent_loss=parent_loss, child_loss=child_loss)


class StructurePrecisionEnum(enum.Enum):
    unknown = -1
    ranging = 1
    exact = 2


def decorate_tree(tree, decorate_value):
    '''
    Transform each node in the tree's id value
    into a tuple of `(decorate_value, node.id)`
    to use common ids across identical graphs but
    discriminable between attached copies.

    Parameters
    ----------
    tree: Glycan
    decorate_value: int

    Returns
    -------
    Glycan
    '''
    for node in tuple(tree):
        node.id = (decorate_value, node.id)
    return tree


def decorated_value(tree):
    '''
    Get the decorating value from a tree's root node.id.

    Returns `None` if the tree is undecorated

    Parameters
    ----------
    tree: Glycan

    Returns
    -------
    int or None
    '''
    try:
        return iter(tree.root.id).next()
    except:
        return None


def get_decorated(tree, id):
    '''
    As :meth:`Glycan.get`, but with awareness of decorated
    node.id attributes. Will fall back to the normal get method
    if the tree is undecorated.

    Parameters
    ----------
    tree: Glycan
    id: int

    Returns
    -------
    Monosaccharide
    '''
    d = decorated_value(tree)
    if d is None:
        return tree.get(id)
    else:
        return tree.get((d, id))


def undecorate_tree(tree):
    '''
    Remove decoration from a tree and reindex its nodes.

    Depends upon  :attr:`Glycan.index` to find nodes

    Parameters
    ----------
    tree: Glycan

    Returns
    -------
    Glycan
    '''
    i = uuid4().int
    for node in tree.index:
        node.id = i
        i += 1
    tree.reindex()
    return tree


class RepeatRecord(object):
    def __init__(self, graph_index, repeat_index, internal_linkage=None,
                 external_linkage=None, multitude=(-1, -1), graph=None,
                 parser=None):
        if graph is None:
            graph = {}
        self.graph_index = graph_index
        self.repeat_index = repeat_index
        self.internal_linkage = internal_linkage
        self.external_linkage = external_linkage
        self.multitude = multitude
        self.graph = {}
        self.original_graph = None
        self.parser = parser

    def __contains__(self, k):
        if self.original_graph is not None:
            return k in self.original_graph
        else:
            return k in self.graph

    def is_exact(self):
        if -1 in self.multitude:
            return StructurePrecisionEnum.unknown
        elif self.multitude[0] == self.multitude[1]:
            return StructurePrecisionEnum.exact
        return StructurePrecisionEnum.ranging

    def expand_inner(self, n=None):
        if n is None:
            if self.multitude[1] != -1:
                n = self.multitude[1]
            elif self.multitude[0] != -1:
                n = self.multitude[0]
            else:
                n = 1

        if self.is_exact() is not StructurePrecisionEnum.unknown:
            if not (self.multitude[0] <= n <= self.multitude[1]):
                raise ValueError("{} is not within the range of {}".format(n, self.multitude))
        sub_unit_indices = sorted(self.graph.keys())

        if self.original_graph is None:
            glycan_graph = Glycan(self.graph[sub_unit_indices[0]], index_method=None).clone()
            self.original_graph = {}
            for k, v in self.graph.items():
                try:
                    self.orignal_graph[k] = glycan_graph.get(v.id)
                except:
                    self.original_graph[k] = self.graph[k]

        else:
            self.graph = {k: v.clone() for k, v in self.graph.items()}
        glycan_graph = Glycan(self.graph[sub_unit_indices[0]], index_method=None)

        graph = {1: glycan_graph.clone(index_method=None)}
        decorate_tree(graph[1], 1)
        parent_residue_index = int(self.internal_linkage["parent_residue_index"])
        parent_atom_replaced = link_replacement_composition_map[self.internal_linkage["parent_atom_replaced"]]
        parent_attachment_position = self.internal_linkage["parent_attachment_position"]
        child_residue_index = int(self.internal_linkage["child_residue_index"])
        child_atom_replaced = link_replacement_composition_map[self.internal_linkage["child_atom_replaced"]]
        child_attachment_position = self.internal_linkage["child_attachment_position"]

        op_stack = []

        for i in range(2, n + 1):
            graph[i] = glycan_graph.clone(index_method=None)
            graph[i] = decorate_tree(graph[i], i)
            parent_graph = graph[i-1]
            child_graph = graph[i]
            parent_node = parent_graph.get((i - 1, parent_residue_index))
            child_node = child_graph.get((i, child_residue_index))
            op_stack.append(
                (form_link, [parent_node, child_node],
                 dict(parent_position=parent_attachment_position,
                      child_position=child_attachment_position,
                      parent_loss=parent_atom_replaced, child_loss=child_atom_replaced)))
        for op in op_stack:
            f, args, kwargs = op
            f(*args, **kwargs)

        self.graph = graph

    def handle_incoming_link(self, parent, parent_position, parent_loss,
                             child_position, child_loss, id=None):
        sub_unit_indices = sorted(self.graph.keys())
        child_graph = self.graph[sub_unit_indices[0]]
        parent = parent()
        child = child_graph.get((1, int(self.external_linkage['child_residue_index'])))
        if parent_loss == Composition("H"):
            child_loss = Composition("OH")

        form_link(parent, child, parent_position=parent_position, child_position=child_position,
                  parent_loss=parent_loss, child_loss=child_loss, id=id)

    def handle_outgoing_link(self, child, parent_position, parent_loss,
                             child_position, child_loss, id=None):
        sub_unit_indices = sorted(self.graph.keys())
        parent_graph = self.graph[sub_unit_indices[-1]]
        child = child()
        parent = parent_graph.get((len(sub_unit_indices), int(self.external_linkage['parent_residue_index'])))
        if isinstance(child, RepeatRecord):
            child.handle_incoming_link(
                lambda: parent, parent_position=parent_position,
                child_position=child_position, parent_loss=parent_loss, child_loss=child_loss, id=id)
        else:
            form_link(
                parent, child, parent_position=parent_position, child_position=child_position,
                parent_loss=parent_loss, child_loss=child_loss, id=id)

    def get_node(self, id, direction=None):
        id = int(id)
        if direction is None or direction == "in":
            sub_unit_indices = sorted(self.graph.keys())
            child_graph = self.graph[sub_unit_indices[0]]
            child = child_graph.get((1, id))
            return child
        elif direction == "out":
            sub_unit_indices = sorted(self.graph.keys())
            parent_graph = self.graph[sub_unit_indices[-1]]
            parent = parent_graph.get((len(sub_unit_indices), id))
            return parent
        else:
            raise Exception("Unknown direction %s" % direction)

    def __root__(self):
        root_node = rootp(self.graph[1])
        if root_node.node_type is Substituent.node_type:
            root_node = root_node.links[1][0].parent
        return root_node

    def prepare_glycan(self):
        glycan = self.graph[1]
        glycan.deindex()
        return glycan


def try_int(v):
    try:
        return int(v)
    except:
        return None


class GlycoCTError(Exception):
    pass


class GlycoCTSectionUnsupported(GlycoCTError):
    pass


class GlycoCT(object):
    '''
    Simple State-Machine parser for condensed GlycoCT representations. Yields
    |Glycan| instances.
    '''

    @classmethod
    def loads(cls, glycoct_str, structure_class=Glycan):
        '''Parse results from |str|'''
        rep = StringIO(glycoct_str)
        return cls(rep, structure_class=structure_class)

    def __init__(self, stream, structure_class=Glycan):
        '''
        Creates a parser of condensed GlycoCT.

        Parameters
        ----------
        stream: basestring or file-like
            A path to a file or a file-like object to be processed
        '''
        self.graph = {}
        self.handle = opener(stream, "r")
        self.state = START
        self.current_repeat = None
        self.in_repeat = False
        self.repeats = {}
        self.postponed = []
        self.root = None
        self._iter = None
        self.structure_class = structure_class

    def _read(self):
        for line in self.handle:
            for token in re.split(r"\s|;", line):
                # logger.debug(token)
                if "" == token.strip():
                    continue
                yield token

    def _reset(self):
        self.graph = {}
        self.root = None
        self.postponed = []

    def __iter__(self):
        '''
        Calls :meth:`parse` and stores it for reuse with :meth:`__next__`
        '''
        self._iter = self.parse()
        return self._iter

    def next(self):
        '''
        Calls :meth:`parse` if the internal iterator has not been instantiated

        '''
        if self._iter is None:
            iter(self)
        return self._iter.next()

    #: Alias for next. Supports Py3 Iterator interface
    __next__ = next

    def deferred_retrieval(self, id, direction=None):
        def retrieval():
            try:
                return self.graph[id]
            except KeyError:
                for repeat_ix, repeater in self.repeats.items():
                    if id in repeater:
                        return repeater.get_node(id, direction)
        return retrieval

    def handle_residue_line(self, line):
        '''
        Handle a base line, creates an instance of |Monosaccharide|
        and adds it to :attr:`graph` at the given index.

        Called by :meth:`parse`
        '''
        _, ix, residue_str = re.split(r"^(\d+)b", line, maxsplit=1)
        residue_dict = res_pattern.search(residue_str).groupdict()

        mods = residue_dict.pop("modifications")
        modifications = OrderedMultiMap()
        if mods is not None:
            for p, mod in modification_pattern.findall(mods):
                modifications[try_int(p)] = modification_map[mod]

        residue_dict["modifications"] = modifications
        is_reduced = "aldi" in modifications[1]
        if is_reduced:
            modifications.pop(1, "aldi")
            residue_dict['reduced'] = True

        conf_stem = residue_dict.pop("conf_stem")
        if conf_stem is not None:
            config, stem = zip(*conf_stem_pattern.findall(conf_stem))
        else:
            config = ('x',)
            stem = ('x',)
        residue_dict['stem'] = stem
        residue_dict['configuration'] = config

        residue_dict["ring_start"], residue_dict["ring_end"] = list(map(
            try_int, residue_dict.pop("indices").split(":")))

        residue_dict['anomer'] = anomer_map[residue_dict['anomer']]
        residue_dict['superclass'] = superclass_map[residue_dict['superclass']]
        residue = monosaccharide.Monosaccharide(**residue_dict)
        if self.in_repeat:
            graph = self.current_repeat.graph
        else:
            graph = self.graph
        graph[ix] = residue

        residue.id = int(ix)
        if self.root is None:
            if self.in_repeat:
                self.root = self.current_repeat
            else:
                self.root = residue

    def handle_residue_substituent(self, line):
        '''
        Handle a substituent line, creates an instance of |Substituent|
        and adds it to :attr:`graph` at the given index. The |Substituent| object is not yet linked
        to a |Monosaccharide| instance.

        Called by :meth:`parse`

        '''
        _, ix, subsituent_str = re.split(r"^(\d+)s:", line, maxsplit=1)
        sub = Substituent(subsituent_str.strip())

        if self.in_repeat:
            graph = self.current_repeat.graph
        else:
            graph = self.graph

        graph[ix] = sub

    def get_node(self, id):
        if self.in_repeat:
            try:
                return self.current_repeat.graph[id]
            except KeyError:
                return self.graph[id]
        else:
            return self.graph[id]

    def handle_linkage(self, line):
        '''
        Handle a linkage line, creates an instance of |Link| and
        attaches it to the two referenced nodes in :attr:`graph`. The parent node is always
        an instance of |Monosaccharide|, and the child node
        may either be an instance of |Monosaccharide| or
        |Substituent| or |Monosaccharide|.

        Called by :meth:`parse`

        See also |Link| for more information on the impact of instantiating
        a |Link| object.
        '''
        id, parent_residue_index, parent_atom_replaced, parent_attachment_position,\
            child_residue_index, child_atom_replaced, child_attachment_position = parse_link(line)

        parent = self.get_node(parent_residue_index)
        child = self.get_node(child_residue_index)

        if isinstance(parent, RepeatRecord):
            self.postponed.append((
                parent.handle_outgoing_link,
                self.deferred_retrieval(child_residue_index),
                parent_attachment_position, parent_atom_replaced,
                child_attachment_position, child_atom_replaced, id))
            return
        if isinstance(child, RepeatRecord):
            self.postponed.append((
                child.handle_incoming_link,
                self.deferred_retrieval(parent_residue_index),
                parent_attachment_position, parent_atom_replaced,
                child_attachment_position, child_atom_replaced, id))
            return

        form_link(
            parent, child,
            parent_position=parent_attachment_position, child_position=child_attachment_position,
            parent_loss=parent_atom_replaced, child_loss=child_atom_replaced, id=id)

    def handle_repeat_stub(self, line):
        match = repeat_line_pattern.search(line).groupdict()
        graph_index = (match['graph_index'])
        repeat_index = (match["repeat_index"])
        repeat = RepeatRecord(graph_index, repeat_index, parser=self)
        self.graph[graph_index] = repeat
        self.repeats[repeat_index] = repeat

    def postprocess(self):
        '''
        Handle all deferred operations such as binding together and expanding
        repeating units. Removes any distinguishing markers on node ids, and
        constructs a new instance of :attr:`structure_class` from the accumulated
        graph

        Returns
        -------
        Glycan
        '''
        for repeat_index, repeater in self.repeats.items():
            repeater.expand_inner()

        for postop in self.postponed:
            logger.debug("Postprocessing %s", postop)
            if debug:
                pdb.set_trace()
            postop[0](*postop[1:])


        return undecorate_tree(self.structure_class(rootp(self.root)))#.reindex()

    def parse(self):
        '''
        Returns an iterator that yields each complete :class:`Glycan` instance
        from the underlying text stream.
        '''
        for line in self._read():
            if RES == line.strip():
                self.state = RES
                # logger.debug("RES")
                if self.root is not None and not self.in_repeat:
                    yield self.postprocess()
                    self._reset()
            elif LIN == line.strip():
                if self.state != RES:
                    raise GlycoCTError("LIN before RES")
                self.state = LIN

            elif REP == line.strip():
                self.state = REP
                logger.debug("REP")
                self.in_repeat = True
                #raise GlycoCTSectionUnsupported(REP)

            elif line.strip()[:3] == REP:
                logger.debug(line)
                if not self.in_repeat:
                    raise GlycoCTError("Encountered {} outside of REP".format(line))
                header_dict = rep_header_pattern.search(line).groupdict()
                repeat_index = (header_dict['repeat_index'])
                repeat_record = self.repeats[repeat_index]
                self.current_repeat = repeat_record
                linkage = internal_link_pattern.search(header_dict['internal_linkage']).groupdict()
                repeat_record.internal_linkage = linkage
                repeat_record.external_linkage = linkage
                repeat_record.multitude = tuple(map(try_int, (header_dict['lower_multitude'],
                                                              header_dict['higher_multitude'])))
                self.state = START
            elif ALT == line.strip():
                raise GlycoCTSectionUnsupported(ALT)
            elif UND == line.strip():
                raise GlycoCTSectionUnsupported(UND)

            elif re.search(r"^(\d+)b", line) and self.state == RES:
                # logger.debug("handling residue")
                self.handle_residue_line(line)

            elif re.search(r"^(\d+)s:", line) and self.state == RES:
                # logger.debug("handling subsituent")
                self.handle_residue_substituent(line)
            elif re.search(r"^(\d+)r:", line) and self.state == RES:
                # raise GlycoCTSectionUnsupported(REP)
                self.handle_repeat_stub(line)
            elif re.search(r"^(\d+):(\d+)", line) and self.state == LIN:
                # logger.debug("handling linkage")
                self.handle_linkage(line)
            else:
                raise GlycoCTError("Unknown format error: {}".format(line))
        self.in_repeat = False
        yield self.postprocess()


def read(stream, structure_class=Glycan):
    '''
    A convenience wrapper for :class:`GlycoCT`
    '''
    return GlycoCT(stream, structure_class=structure_class)


def loads(glycoct_str, structure_class=Glycan):
    '''
    A convenience wrapper for :meth:`GlycoCT.loads`

    As additional convenience, this function does not return an
    iterator over glycans, and returns a single instance if only
    one is present, or a list of instances otherwise.
    '''

    g = GlycoCT.loads(glycoct_str, structure_class=structure_class)
    first = g.next()
    second = None
    try:
        second = g.next()
        collection = [first, second] + list(g)
        return collection
    except StopIteration:
        return first