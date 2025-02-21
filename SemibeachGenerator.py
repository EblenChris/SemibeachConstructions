from ArcDiagrams import *


class SemibeachGenerator:
    """
    Generates a random semibeach.

    Attributes:
        n (int): The number of possible endpoints for arcs, which are 1, 2, ..., n.
        k (float): A parameter controlling the size of the beach.
        q (int): The size of the finite field.
        row_size (int): The number of rows, determined by floor(k).
        col_size (int): The number of columns, determined by 2k + 1.
        fin_field_units (List): The nonzero elements of the finite field.
        fin_field (List): The elements of the finite field including zero.
        arcs (List): A list of possible arcs in the diagram.
    """
    def __init__(self, n: int, k: float, q: int):
        """
        Initializes the SemibeachGenerator with given Args.

        Args:
            n (int): The number of possible endpoints for arcs, which are 1, 2, ..., n.
            k (float): Determines the size of the beach.
            q (int): The order of the finite field.
        """
        self.n = n
        self.k = k
        self.q = q
        self.row_size = floor(k)
        self.col_size = int(2 * k + 1)
        self.fin_field_units = list(symbols('f1:%d' % q)) if q != 2 else [1]
        self.fin_field = [0] + self.fin_field_units
        self.arcs = list(((i, j), a) for i in range(1, n + 1) for j in range(i + 1, n + 1)
                         for a in self.fin_field_units) + list(((i, i), 0) for i in range(1, n + 1))

    def seed_choices(self, arc: Tuple[Tuple[int, int], Any], j: int) -> List[Tuple[Tuple[int, int], Any]]:
        """
        Returns a list of left or right seed choices based on the given arc and column index.

        Args:
            arc: The arc to consider.
            j (int): The column index.

        Returns:
            List: A list of valid seed arcs.
        """
        if j % 2 == 0:
            return list(arc_1 for arc_1 in self.arcs if Arc(arc_1).re == self.n)
        else:
            return list(arc_1 for arc_1 in self.arcs if Arc(arc_1).le == Arc(arc).le and
                        (Arc(arc_1).re != self.n or Arc(arc).le == self.n))

    def conflict_resolutions(self, arc: Tuple[Tuple[int, int], Any], j: int) -> List[Tuple[Tuple[int, int], Any]]:
        """
        Determines potential conflict resolutions for the given arc and column index.

        Args:
            arc: The arc to resolve conflicts for.
            j (int): The column index.

        Returns:
            List: A list of possible conflict resolutions.
        """
        if j % 2 == 0:
            return list(arc_1 for arc_1 in self.arcs if Arc(arc_1).re == Arc(arc).re and Arc(arc_1).le > Arc(arc).le)
        else:
            return list(arc_1 for arc_1 in self.arcs if Arc(arc_1).le == Arc(arc).le and Arc(arc_1).re < Arc(arc).re)

    def random_semibeach(self, want_beach: bool) -> Tuple[List[List[List[Tuple[Tuple[int, int], Any]]]], List[Tuple]]:
        """
        Generates a random semibeach. The resulting semibeach will be a beach if the boolean want_beach is True.

        Args:
            want_beach (bool): If True, forces the output to be a beach; otherwise, it may be a semibeach.

        Returns:
            Tuple: A tuple containing the generated semibeach and its shifts.
        """
        b = [[list() for _ in range(self.col_size)] for _ in range(self.row_size)]
        shifts = []

        # initial left seed choice
        b[0][2].append(random.choice(self.seed_choices((), 2)))
        shifts.append((Wave(b[0][2]).face, 2))

        def conflict_pair(semibeach, col: int):
            for row in range(int((col - 2) // 2) + 1):
                f_1 = Wave(semibeach[row][col]).face
                for h in range(row + 1, int((col - 2) // 2) + 1):
                    f_2 = Wave(semibeach[h][col]).face
                    if (ArcDiagram([f_1, f_2], self.n, self.q).is_wave()
                            and Arc(f_1).dim > -1 and Arc(f_2).dim > -1):
                        if Arc(f_1).arc_less_eq(f_2):
                            return [(h, f_2), (row, f_1)]
                        else:
                            return [(row, f_1), (h, f_2)]
            return None

        def losing_arc(c):
            return c[1][1]

        def winning_arc(c):
            return c[0][1]

        def losing_row(c):
            return c[1][0]

        def winning_row(c):
            return c[0][0]

        for j in range(3, self.col_size):

            for i in range((j - 2) // 2 + 1):
                if len(b[i][j - 1]) > 0:
                    b[i][j].append(Wave(b[i][j - 1]).face)

            # add a seed for swell j
            f = Wave(b[(j - 2) // 2][j - 1]).face if j % 2 == 1 else []
            b[(j - 2) // 2][j].append(random.choice(self.seed_choices(f, j)))
            shifts.append((Wave(b[(j - 2) // 2][j]).face, j)) if j % 2 == 0 else (
                shifts.append(([Wave(b[(j - 2) // 2][j]).face, Wave(b[(j - 2) // 2][j - 1]).face], j)))

            # get the union of the nontrivial faces in swell j
            union_faces = [Wave(b[i][j]).face for i in range((j - 2) // 2 + 1) if Arc(Wave(b[i][j]).face).dim > -1]

            # coin_flip determines if we keep resolving conflicts or not
            coin_flip = 1 if j < self.col_size - 1 or want_beach else random.choice([0, 1])

            while coin_flip == 1 and ArcDiagram(union_faces, self.n, self.q).is_set_partition() is False:

                # randomly pick a resolving arc
                resolving_arc = random.choice(self.conflict_resolutions(losing_arc(conflict_pair(b, j)), j))

                # add the conflict_pair resolution to the shift set
                shifts.append(((winning_arc(conflict_pair(b, j)), winning_row(conflict_pair(b, j))),
                               ([losing_arc(conflict_pair(b, j)), resolving_arc], losing_row(conflict_pair(b, j))), j))

                # update the union of nontrivial faces
                union_faces.remove(losing_arc(conflict_pair(b, j)))
                if Arc(resolving_arc).dim > -1:
                    union_faces.append(resolving_arc)

                # add the resolving arc to the losing row
                b[losing_row(conflict_pair(b, j))][j].append(resolving_arc)

                # re-flip coin if we are in the final swell
                coin_flip = 1 if j < self.col_size - 1 or want_beach else random.choice([0, 1])

        return b, shifts


class BeachGenerator:
    """
    A class for generating random set partitions and random beaches with a fixed final set partition.

    Attributes:
        n (int): The number of possible endpoints for arcs, which are 1, 2, ..., n.
        k (float): Parameter controlling the size of the beach.
        q (int): The order of the finite field.
        row_size (int): The number of rows in the beach.
        col_size (int): The number of columns in the beach.
        fin_field_units (List): List of nonzero elements in the finite field.
        fin_field (List): List of all elements in the finite field, including zero.
        arcs (List): List of all possible arcs in the diagram.
    """
    def __init__(self, n: int, k: float, q: int):
        """
        Initializes the BeachGenerator with given Args.

        Args:
            n (int): The number of possible endpoints for arcs, which are 1, 2, ..., n.
            k (float): Parameter controlling the size of the beach.
            q (int): The order of the finite field.
        """
        self.n = n
        self.k = k
        self.q = q
        self.row_size = floor(k)
        self.col_size = int(2 * k + 1)
        self.fin_field_units = list(symbols('f1:%d' % q)) if q != 2 else [1]
        self.fin_field = [0] + self.fin_field_units
        self.arcs = list(((i, j), a) for i in range(1, n + 1) for j in range(i + 1, n + 1)
                         for a in self.fin_field_units) + list(((i, i), 0) for i in range(1, n + 1))

    @property
    @lru_cache(maxsize=None)  # Cache the result to avoid recalculation
    def set_partitions_of_n(self):
        """
        Retrieves all (unlabeled) set partitions of the set {1, 2, ..., n}.

        Returns:
            tuple: A tuple of set partitions, where each partition is represented as a tuple of sets.
        """
        return tuple(SetPartitionGenerator(self.n, 2).labeled_set_partitions())

    @property
    @lru_cache(maxsize=None)  # Cache the result to avoid recalculation
    def set_partitions_of_n_minus_1(self):
        """
        Retrieves all (unlabeled) set partitions of the set {1, 2, ..., n - 1}.

        Returns:
            tuple: A tuple of set partitions, where each partition is represented as a tuple of sets.
        """
        return tuple(SetPartitionGenerator(self.n - 1, 2).labeled_set_partitions())

    def random_super_character(self):
        """
        Generates a random labeled set partition corresponding to a supercharacter in the kth induction/restriction step.

        Returns:
            List: A randomly generated set partition.
        """
        return random.choice([sp for sp in self.set_partitions_of_n if len(sp) <= floor(self.k)]) \
            if int(2 * self.k) % 2 == 0 \
            else random.choice([sp for sp in self.set_partitions_of_n_minus_1 if len(sp) <= floor(self.k)])

    def random_beach(self, set_partition: List[Tuple[Tuple[int, int], Any]]) \
            -> Tuple[List[List[List[Tuple[Tuple[int, int], Any]]]], List[Tuple]]:
        """
        Generates a random beach for a given labeled set partition.

        Args:
            set_partition (list): The set partition to visualize.

        Returns:
            tuple: A 2D list representing the beach and a list of transformations applied.
        """
        b = [[list() for _ in range(self.col_size)] for _ in range(self.row_size)]
        shifts = []

        def append_arcs(beach_in_progress, arcs_with_rows_in_progress, swell):
            if not arcs_with_rows_in_progress:
                # assign the arcs of the given set partition randomly to the rows in the last swell.
                unused_rows = [i for i in range(self.row_size - 1)] if int(2 * self.k) % 2 == 0 else [i for i in range(
                    self.row_size)]
                for arc in set_partition:
                    if Arc(arc).re != self.n:
                        picked_row = random.choice(unused_rows)
                        arcs_with_rows_in_progress.append((picked_row, arc))
                        unused_rows.remove(picked_row)
                    else:
                        arcs_with_rows_in_progress.append((self.row_size - 1, arc))
                if self.n not in ArcDiagram(set_partition, self.n, self.q).rights and int(2 * self.k) % 2 == 0:
                    arcs_with_rows_in_progress.append((self.row_size - 1, ((self.n, self.n), 0)))

            for i, arc in arcs_with_rows_in_progress:
                beach_in_progress[i][swell].append(arc)

        def moved_arc_choices(beach_in_progress, arcs_with_rows_in_progress, conflict_pair, swell):
            seed_face_row, seed_face = (swell - 2) // 2, Wave(beach_in_progress[(swell - 2) // 2][swell]).face
            if conflict_pair is None:
                curr_set_partition = [arc for i, arc in arcs_with_rows_in_progress if Arc(arc).dim > -1]
                if swell % 2 == 0:
                    dead_arc_choices = [(i, ((right, right), 0)) for i in range((swell - 2) // 2) for right in
                                        range(Arc(seed_face).le + 1, Arc(seed_face).re) if
                                        len(beach_in_progress[i][swell]) == 0 and right not in
                                        ArcDiagram(curr_set_partition, self.n, self.q).rights]
                else:
                    dead_arc_choices = [(i, ((left, left), 0)) for i in range((swell - 2) // 2) for left in
                                        range(Arc(seed_face).le + 1, Arc(seed_face).re) if
                                        len(beach_in_progress[i][swell]) == 0 and left not in
                                        ArcDiagram(curr_set_partition, self.n, self.q).lefts]
                if Arc(seed_face).dim > -1:
                    choices = ([(i, arc) for i, arc in arcs_with_rows_in_progress if
                                Arc(seed_face).le <= Arc(arc).le < Arc(arc).re <= Arc(seed_face).re] + dead_arc_choices)
                else:
                    choices = [((swell - 2) // 2, seed_face)]
            else:
                just_won_arc_row, just_won_arc = conflict_pair[0]
                just_moved_arc_row, just_moved_arc = conflict_pair[1]
                if swell % 2 == 0:
                    choices = [(seed_face_row, seed_face)] if seed_face == just_won_arc else conflict_pair
                else:
                    if Arc(just_won_arc).le == Arc(seed_face).le:
                        choices = [(just_won_arc_row, just_won_arc)]
                    elif Arc(just_moved_arc).le == Arc(seed_face).le:
                        choices = [(just_moved_arc_row, just_moved_arc)]
                    else:
                        choices = [(i, arc) for i, arc in conflict_pair if Arc(seed_face).le < Arc(arc).le]
            return choices

        def winning_arc_choices(beach_in_progress, arcs_with_rows_in_progress, moved_arc, swell):
            seed_face = Wave(beach_in_progress[(swell - 2) // 2][swell]).face

            if swell % 2 == 0:
                choices = [(i, arc) for i, arc in arcs_with_rows_in_progress if
                           Arc(seed_face).le <= Arc(arc).le < Arc(moved_arc).le and Arc(moved_arc).re < Arc(arc).re]
            else:
                if Arc(moved_arc).le == Arc(seed_face).le:
                    choices = [(i, arc) for i, arc in arcs_with_rows_in_progress if Arc(arc).le < Arc(seed_face).le
                               and Arc(moved_arc).re < Arc(arc).re]
                else:
                    choices = [(i, arc) for i, arc in arcs_with_rows_in_progress if Arc(arc).le < Arc(moved_arc).le
                               and Arc(moved_arc).re < Arc(arc).re <= Arc(seed_face).re]

            return choices

        # Given the faces in swell j, start moving backward to the seed for swell j
        def resolve_swell(arcs_with_rows, conflict_pair, j):

            if not conflict_pair:
                # place the given arcs in their respective rows in swell j
                append_arcs(b, arcs_with_rows, j)

                if j % 2 == 1 and not b[(j - 2) // 2][j]:
                    curr_set_partition = [arc for i, arc in arcs_with_rows if Arc(arc).dim > -1]
                    seed_face = random.choice([((left, left), 0) for left in range(1, self.n + 1) if left not in
                                               ArcDiagram(curr_set_partition, self.n, self.q).lefts])
                    for h in range(j, self.col_size):
                        b[(j - 2) // 2][h].append(seed_face)
                    arcs_with_rows.append(((j - 2) // 2, seed_face))

            moved_arc_row, moved_arc = random.choice(moved_arc_choices(b, arcs_with_rows, conflict_pair, j))

            if (moved_arc_row, moved_arc) not in arcs_with_rows:
                for h in range(j, self.col_size):
                    b[moved_arc_row][h].append(moved_arc)
                arcs_with_rows.append((moved_arc_row, moved_arc))

            if j % 2 == 0:
                if moved_arc_row == (j - 2) // 2:
                    arcs_with_rows.remove(((j - 2) // 2, moved_arc))
                    shifts.append((moved_arc, j))
                    if j > 2:
                        resolve_swell(arcs_with_rows, None, j - 1)
                else:
                    winning_arc_row, winning_arc = random.choice(winning_arc_choices(b, arcs_with_rows, moved_arc, j))
                    origin_arc = ((Arc(winning_arc).le, Arc(moved_arc).re), random.choice(self.fin_field_units))
                    b[moved_arc_row][j].append(origin_arc)
                    arcs_with_rows.remove((moved_arc_row, moved_arc))
                    arcs_with_rows.append((moved_arc_row, origin_arc))
                    shifts.append(((winning_arc, winning_arc_row), ([origin_arc, moved_arc], moved_arc_row), j))
                    resolve_swell(arcs_with_rows, [(winning_arc_row, winning_arc), (moved_arc_row, origin_arc)], j)
            elif j % 2 == 1 and Arc(moved_arc).le < self.n:
                winning_arc_row, winning_arc = random.choice(winning_arc_choices(b, arcs_with_rows, moved_arc, j)) \
                    if winning_arc_choices(b, arcs_with_rows, moved_arc, j) else (None, None)
                if winning_arc is None:
                    origin_arc = [(Arc(moved_arc).le, self.n), random.choice(self.fin_field_units)]
                    b[moved_arc_row][j].append(origin_arc)
                    arcs_with_rows.remove((moved_arc_row, moved_arc))
                    arcs_with_rows.append((moved_arc_row, origin_arc))
                    shifts.append(([origin_arc, moved_arc], j))
                    resolve_swell(arcs_with_rows, None, j - 1)
                else:
                    origin_arc = [(Arc(moved_arc).le, Arc(winning_arc).re), random.choice(self.fin_field_units)]
                    b[moved_arc_row][j].append(origin_arc)
                    arcs_with_rows.remove((moved_arc_row, moved_arc))
                    arcs_with_rows.append((moved_arc_row, origin_arc))
                    shifts.append(((winning_arc, winning_arc_row), ([origin_arc, moved_arc], moved_arc_row), j))
                    resolve_swell(arcs_with_rows, [(winning_arc_row, winning_arc), (moved_arc_row, origin_arc)], j)
            elif j % 2 == 1 and Arc(moved_arc).le == self.n:
                shifts.append(([moved_arc, moved_arc], j))
                resolve_swell(arcs_with_rows, None, j - 1)

        resolve_swell([], None, self.col_size - 1)

        return b, shifts[::-1]
