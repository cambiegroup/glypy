import json
from collections import namedtuple, defaultdict, deque

from glypy.io import glycoct
from glypy.structure.glycan_composition import HashableGlycanComposition


EnzymeEdge = namedtuple("EnzymeEdge", ("parent", "child", "enzyme"))


def _enzyme_graph_inner():
    return defaultdict(set)


class EnzymeGraph(object):
    def __init__(self, graph, seeds=None):
        self.graph = graph
        self.seeds = set()
        if seeds is None:
            seeds = self.parentless()
        self.seeds.update(seeds)

    def clone(self):
        graph = defaultdict(_enzyme_graph_inner)
        for outer_key, outer_value in self.graph.items():
            for inner_key, inner_value in outer_value.items():
                graph[outer_key][inner_key] = inner_value.copy()
        return self.__class__(graph, self.seeds.copy())

    def nodes(self):
        acc = set()
        acc.update(self.graph)
        for i, v in enumerate(self.graph.values()):
            acc.update(v)
        return acc

    def edges(self):
        edges = set()
        for outer_key, outer_value in self.graph.items():
            for inner_key, inner_value in outer_value.items():
                for val in inner_value:
                    edges.add(EnzymeEdge(outer_key, inner_key, val))
        return edges

    def node_count(self):
        acc = set()
        acc.update(self.graph)
        for i, v in enumerate(self.graph.values()):
            acc.update(v)
        return len(acc)

    def edge_count(self):
        edges = 0
        for outer_key, outer_value in self.graph.items():
            for inner_key, inner_value in outer_value.items():
                edges += len(inner_value)
        return edges

    def __repr__(self):
        return "{}({:d})".format(self.__class__.__name__, self.node_count())

    def enzymes(self):
        enzyme_set = set()
        for outer_key, outer_value in self.graph.items():
            for inner_key, inner_value in outer_value.items():
                enzyme_set.update(inner_value)
        return enzyme_set

    def remove_enzyme(self, enzyme):
        edges_removed = list()
        for outer_key, outer_value in list(self.graph.items()):
            for inner_key, inner_value in list(outer_value.items()):
                if enzyme in inner_value:
                    inner_value.remove(enzyme)
                    edges_removed.append((outer_key, inner_key))
                if not inner_value:
                    outer_value.pop(inner_key)
                if not outer_value:
                    self.graph.pop(outer_key)
        nodes_to_remove = self.parentless() - self.seeds
        while nodes_to_remove:
            for node in nodes_to_remove:
                self.remove(node)
            nodes_to_remove = self.parentless() - self.seeds
        return edges_removed

    def parents(self, target):
        parents = []
        for outer_key, outer_value in self.graph.items():
            for inner_key, inner_value in outer_value.items():
                if inner_key == target:
                    parents.append(outer_key)
        return parents

    def parentless(self):
        is_parent = set(self.graph)
        is_parented = set()
        for i, v in enumerate(self.graph.values()):
            is_parented.update(v)
        return is_parent - is_parented

    def children(self, target):
        children = []
        children.extend(self.graph[target])
        return children

    def remove(self, prune):
        items = deque([prune])
        i = 0
        while items:
            node = items.popleft()
            if node in self.graph:
                i += 1
                self.graph.pop(node)
        return i

    def _dump(self):
        data_structure = {
            "seeds": [str(sd) for sd in self.seeds],
            "enzymes": list(self.enzymes()),
            "graph": {}
        }
        outgraph = {}
        for outer_key, outer_value in self.graph.items():
            outgraph_inner = dict()
            for inner_key, inner_value in outer_value.items():
                outgraph_inner[str(inner_key)] = list(inner_value)
            outgraph[str(outer_key)] = outgraph_inner
        data_structure['graph'] = outgraph
        return data_structure

    def dump(self, fh):
        d = self._dump()
        json.dump(d, fh, sort_keys=True, indent=2)

    def dumps(self):
        d = self._dump()
        return json.dumps(d, sort_keys=True, indent=2)

    @classmethod
    def _load_entity(self, entity):
        return entity

    @classmethod
    def _load(cls, data_structure):
        seeds = {cls._load_entity(sd) for sd in data_structure["seeds"]}
        graph = dict()
        for outer_key, outer_value in data_structure["graph"].items():
            outgraph_inner = dict()
            for inner_key, inner_value in outer_value.items():
                outgraph_inner[cls._load_entity(inner_key)] = set(inner_value)
            graph[cls._load_entity(outer_key)] = outgraph_inner
        inst = cls(graph, seeds)
        return inst

    @classmethod
    def loads(cls, text):
        data = json.loads(text)
        return cls._load(data)

    @classmethod
    def load(cls, fd):
        data = json.load(fd)
        return cls._load(data)

    def __eq__(self, other):
        return self.graph == other.graph

    def __ne__(self, other):
        return self.graph != other.graph


# This may be too memory intensive to use on large graphs because
# a single :class:`~.Glycan` instance uses many times the memory that
# a :class:`~.GlycanComposition` does.
class GlycanStructureEnzymeGraph(EnzymeGraph):

    @classmethod
    def _load_entity(self, entity):
        return glycoct.loads(entity)


class GlycanCompositionEnzymeGraph(EnzymeGraph):

    @classmethod
    def _load_entity(self, entity):
        return HashableGlycanComposition.parse(entity)
