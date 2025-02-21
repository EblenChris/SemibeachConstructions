from sympy import *
from math import comb
from typing import *
from functools import lru_cache
import random


class Arc:
    """
    Represents a single arc within an arc diagram. Contains methods to extract endpoints, label, and other arc-specific
    properties.

    Attributes:
        arc (Tuple[Tuple[int, int], Any]): A tuple containing the arc's endpoints and its label.
        le (int): The left endpoint of the arc.
        re (int): The right endpoint of the arc.
        lab (Any): The label of the arc.
        dim (int): The dimension or length of the arc, defined as the difference between the right and left endpoints.

    Methods:
        arc_less_eq(arc_2) -> bool:
            Checks whether an arc is less than or equal to arc_2 in terms of dimension.
    """
    def __init__(self, arc: Tuple[Tuple[int, int], Any]):
        """
        Initializes an Arc instance.

        Args:
            arc (Tuple[Tuple[int, int], Any]): A tuple containing the endpoints and label of the arc.
        """
        self.arc = arc
        self.le = arc[0][0]
        self.re = arc[0][1]
        self.lab = arc[1]
        self.dim = self.re - self.le - 1

    def __repr__(self):
        return f'Arc(arc={self.arc})'

    def arc_less_eq(self, arc_2: Tuple[Tuple[int, int], Any]) -> bool:
        """
        Compares two arcs based on their dimension (length).

        Args:
            arc_2 (Tuple[Tuple[int, int], Any]): The second arc to compare against.

        Returns:
            bool: True if the current arc is less than or equal to the second arc in dimension, False otherwise.
        """
        return self.dim <= Arc(arc_2).dim


class Wave:
    """
    Represents a wave in an arc diagram, defined as a collection of arcs with certain properties.
    The wave includes the set of left endpoints, right endpoints, and the smallest and largest arcs.

    Attributes:
        wave (List[Tuple[Tuple[int, int], Any]]): A list of arcs in the wave, where each arc is represented
            as a tuple containing the endpoints and a label.
        face (Tuple[Tuple[int, int], Any]): The smallest arc in the wave, based on its dimension (distance
            between endpoints). If the wave is empty, `face` is an empty list.
        back (Tuple[Tuple[int, int], Any]): The largest arc in the wave, based on its dimension. If the wave
            is empty, `back` is an empty list.
        lefts (Set[int]): A set of all left endpoints of arcs in the wave.
        rights (Set[int]): A set of all right endpoints of arcs in the wave.

    Methods:
        is_left() -> bool:
            Checks if a wave is a left wave, i.e. has at most one left endpoint.

        is_right() -> bool:
            Checks if a wave is a right wave, i.e. has at most one right endpoint.
    """
    def __init__(self, wave: List[Tuple[Tuple[int, int], Any]]):
        """
        Initializes a Wave instance.

        Args:
            wave (List[Tuple[Tuple[int, int], Any]]): A list of arcs in the wave, where each arc is represented
                as a tuple containing the endpoints and a label.
        """
        self.wave = sorted(wave, key=lambda sub: abs(sub[0][1] - sub[0][0]))
        self.face = self.wave[0] if len(self.wave) > 0 else []
        self.back = self.wave[-1] if len(self.wave) > 0 else []
        self.lefts = {Arc(arc).le for arc in wave}
        self.rights = {Arc(arc).re for arc in wave}

    def __repr__(self):
        return f'Wave(wave={self.wave})'

    def is_left(self) -> bool:
        """
        Checks if the wave is a left wave. A left wave has at most one unique left endpoint.

        Returns:
            bool: True if the wave is a left wave, False otherwise.
        """
        return len(self.lefts) <= 1

    def is_right(self) -> bool:
        """
        Checks if the wave is a right wave. A right wave has at most one unique right endpoint.

        Returns:
            bool: True if the wave is a right wave, False otherwise.
        """
        return len(self.rights) <= 1


class SetPartitionGenerator:
    """
    A class to handle set partitions, labeled set partitions, and the computation
    of the generalized Bell number for given values of n and q.

    Attributes:
        n (int): The number of possible endpoints for arcs, which are 1, 2, ..., n.
        q (int): The size of the finite field used for labeling set partitions.

    Methods:
        bell() -> int:
           Computes the generalized Bell number for the given n and q, which represents the number
           of labeled set partitions.

        labeled_set_partitions() -> Generator[List[Tuple[Tuple[int, int], Any]], None, None]:
           Yields the labeled set partitions of {1, 2, ..., n}.
    """
    def __init__(self, n: int, q: int):
        """
        Initializes the SetPartitionGenerator with the set size and finite field size.

        Args:
            n (int): The size of the set {1, 2, ..., n}.
            q (int): The size of the finite field used for labeling partitions.
        """
        self.n = n
        self.q = q
        self.fin_field_units = list(symbols('f1:%d' % q)) if q != 2 else [1]
        self.fin_field = [0] + self.fin_field_units

    def __repr__(self):
        return f'SetPartitionGenerator(n={self.n}, q={self.q})'

    def bell(self) -> int:
        """
        Computes the generalized Bell number for the given n and q, which represents the number of labeled set
        partitions. This method uses a recursive approach with caching to optimize performance for repeated subproblem
        calculations.

        Returns:
            int: The generalized Bell number.
        """

        @lru_cache(None)
        def recursion_bell(m):
            if m == 0:
                return 1
            return sum(comb(m - 1, i) * (self.q - 1) ** (m - 1 - i) * recursion_bell(i) for i in range(m))

        return recursion_bell(self.n)

    def labeled_set_partitions(self) -> Generator[List[Tuple[Tuple[int, int], Any]], None, None]:
        """
        Generates labeled set partitions of {1, 2, ..., n}. Each set partition is translated into a list of arcs, where
        each arc is labeled with an element of the finite field.  This implementation uses a generator to avoid memory
        overhead.

        Returns:
            Generator[List[Tuple[Tuple[int, int], Any]]]: A generator that yields labeled set partitions of
            {1, 2, ..., n} one at a time.
        """
        @lru_cache(None)
        def make_set_partitions(m):
            if m == 0:
                yield []
            else:
                for p in make_set_partitions(m - 1):
                    yield [[m]] + p
                    for i, s in enumerate(p):
                        yield p[:i] + [s + [m]] + p[i + 1:]

        set_parts = make_set_partitions(self.n)
        for set_part in set_parts:
            diagram_set_part = []
            for part in set_part:
                if len(part) > 1:
                    for i in range(len(part) - 1):
                        arc = (part[i], part[i + 1])
                        label = random.choice(list(self.fin_field_units))
                        diagram_set_part.append((arc, label))
            yield diagram_set_part


class ArcDiagram:
    """
   Represents an arc diagram with arcs labeled by a finite field. Supports various
   operations such as checking structural properties, generating labeled set partitions,
   and computing related combinatorial numbers.

   Attributes:
       arc_diagram (List[Tuple[Tuple[int, int], Any]]): A list of arcs, where each arc is a tuple
           consisting of a pair of integers (the endpoints) and a label.
       n (int): The number of possible endpoints for arcs, which are 1, 2, ..., n.
       q (int): The order of the finite field for labeling the arcs.
       lefts (set[int]): The set of all left endpoints of arcs in the diagram.
       rights (set[int]): The set of all right endpoints of arcs in the diagram.

   Methods:
       set_partitions():
           Generates all set partitions of the set {1, 2, ..., n}.

       labeled_set_partitions():
           Converts the set partitions of {1, 2, ..., n} into labeled set partitions.

       is_wave() -> bool:
           Checks whether the arc diagram is a wave (one of the sets of left or right endpoints has at most one
           element).

       is_set_partition() -> bool:
           Checks whether the arc diagram corresponds to a valid set partition.

       is_cascade() -> bool:
           Checks whether the arc diagram is a cascade of waves.

       crossings(arc_diagram_2) -> set:
           Finds the crossing pairs of arcs between this diagram and another arc diagram.

       nestings(arc_diagram_2) -> set:
           Finds the nesting pairs of arcs between this diagram and another arc diagram.

       bell() -> int:
           Computes the generalized Bell number for the given n and q, which represents the number
           of labeled set partitions.
       """
    def __init__(self, arc_diagram: List[Tuple[Tuple[int, int], Any]], n: int, q: int):
        """
        Initializes an ArcDiagram instance.

        Args:
            arc_diagram (List[Tuple[Tuple[int, int], Any]]): List of arcs, each represented by
                a tuple containing the endpoints and the label.
            n (int): The number of possible endpoints for arcs, which are 1, 2, ..., n.
            q (int): The order of the finite field for labeling arcs.
        """
        self.arc_diagram = arc_diagram
        self.n = n
        self.q = q
        self.lefts = {Arc(arc).le for arc in arc_diagram}
        self.rights = {Arc(arc).re for arc in arc_diagram}

    def __repr__(self):
        return f'ArcDiagram(n={self.n}, q={self.q}, arcs={self.arc_diagram})'

    def is_wave(self) -> bool:
        """
        Checks whether the arc diagram is a wave. A wave is defined as having at least one of the sets of left or right
        endpoints with at most one element.

        Returns:
            bool: True if the arc diagram is a wave, False otherwise.
        """
        return len(self.lefts) <= 1 or len(self.rights) <= 1

    def is_set_partition(self) -> bool:
        """
        Checks whether the arc diagram is a labeled set partition. A labeled set partition is an arc diagram with no
        trivial arcs (arcs where the left and right endpoint are the same) in which no two arcs form a wave.

        Returns:
            bool: True if the arc diagram is a labeled set partition, False otherwise.
        """
        return all(arc_1 == arc_2 or not ArcDiagram([arc_1, arc_2], self.n, self.q).is_wave() for arc_1 in
                   self.arc_diagram for arc_2 in self.arc_diagram)

    def is_cascade(self) -> bool:
        """
        Checks whether the arc diagram is a cascade of waves. A cascade of waves is defined by the arcs' order and the
        alternating wave properties.

        Returns:
            bool: True if the arc diagram is a cascade of waves, False otherwise.
        """
        arc_diagram = sorted(self.arc_diagram, key=lambda sub: abs(sub[0][1] - sub[0][0]), reverse=True)
        if len(arc_diagram) != 0 and Arc(arc_diagram[0]).re != self.n:
            return False
        else:
            for i in range(1, len(arc_diagram)):
                if any([not ArcDiagram([arc_diagram[i], arc_diagram[i - 1]], self.n, self.q).is_wave(),
                       i % 2 == 1 and Wave([arc_diagram[i], arc_diagram[i - 1]]).is_right(),
                       i % 2 == 0 and Wave([arc_diagram[i], arc_diagram[i - 1]]).is_left()]):
                    return False
            return True

    def crossings(self, arc_diagram_2: List[Tuple[Tuple[int, int], Any]]) -> Set:
        """
        Finds the crossing pairs of arcs between this diagram and another arc diagram. A crossing occurs if one arc's
        left endpoint is between the left and right endpoints of another arc and its right endpoint is to the right of
        both.

        Args:
            arc_diagram_2 (ArcDiagram): The second arc diagram to check for crossings.

        Returns:
            Set: A set of tuple pairs, where each pair consists of two arcs that cross. The first arc in each tuple is
             the arc whose left endpoint is more left.
        """
        return set((arc_1, arc_2) for arc_1 in self.arc_diagram for arc_2 in arc_diagram_2
                   if Arc(arc_1).le < Arc(arc_2).le < Arc(arc_1).re < Arc(arc_2).re)

    def nestings(self, arc_diagram_2: List[Tuple[Tuple[int, int], Any]]) -> Set:
        """
        Finds the nesting pairs of arcs between this diagram and another arc diagram. A nesting occurs if one arc is
        fully contained within another, based on their endpoints.

        Args:
            arc_diagram_2 (ArcDiagram): The second arc diagram to check for nestings.

        Returns:
            Set: A set of tuple pairs, where each pair consists of two arcs that nest within each other. The first arc
            in each tuple is the larger arc.
        """
        return set((arc_1, arc_2) for arc_1 in self.arc_diagram for arc_2 in arc_diagram_2
                   if Arc(arc_1).le < Arc(arc_2).le < Arc(arc_2).re < Arc(arc_1).re)

