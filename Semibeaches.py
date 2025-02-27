from SemibeachGenerator import *


class Semibeach:
    """
    Represents a semibeach structure, which is a combinatorial object related to arc diagrams and partitions.
    """
    def __init__(self, semibeach: Tuple[List[List[List[Tuple[Tuple[int, int], Any]]]], List[Tuple]], n, k, q: int):
        """
        Initializes a Semibeach object.

        Args:
            semibeach (Tuple): A tuple containing a matrix representation of the semibeach and associated shifts.
            q (int): The order of the finite field.
        """
        self.semibeach = semibeach[0]
        self.shifts = semibeach[1]
        self.n = n
        self.k = k
        self.q = q
        self.path = [[Wave(self.semibeach[i][j]).face for i in range((j - 2) // 2 + 1) if len(self.semibeach[i][j]) > 0
                      and Arc(Wave(self.semibeach[i][j]).face).dim > -1] for j in range(len(self.semibeach[0]))]

    @property
    def semibeach_labels(self) -> List[List[List[Any]]]:
        """
        Computes labels for the semibeach.

        Returns:
            List: An array whose entries contains the labels associated with the corresponding entry of the semibeach.
            The list of labels in entry (i, j) corresponds to the labels of the arcs in entry (i, j) of the semibeach
            in descending arc order.
        """
        labels = [[[] for _ in range(int(2 * self.k) + 1)] for _ in range(floor(self.k))]
        for i in range(floor(self.k)):
            for j in range(int(2 * self.k) + 1):
                if j == 2 * i + 3 and ((self.n, self.n), 0) in self.semibeach[i][j]:
                    labels[i][j].append(0)
                else:
                    for h in range(len(self.semibeach[i][j])):
                        labels[i][j].append(self.semibeach[i][j][h][1])
                labels[i][j].reverse()
        return labels

    def cascade_and_partition_checker(self) -> bool:
        """
        Checks whether the structure satisfies the cascade and set partition conditions. Prints diagnostic messages
        indicating where the conditions fail.

        Returns:
            bool: True if every two successive swells" symmetric difference is a cascade of waves and every swell"s
            union of nontrivial faces is a set partition, False otherwise.
        """
        for j in range(1, int(2 * self.k) + 1):
            if (not ArcDiagram([arc for arc in self.path[j] if arc not in self.path[j - 1]] +
                               [arc for arc in self.path[j - 1] if arc not in self.path[j]], self.n, self.q)
                    .is_cascade() or not ArcDiagram(self.path[j], self.n, self.q).is_set_partition()):
                if not (ArcDiagram([arc for arc in self.path[j] if arc not in self.path[j - 1]] +
                                   [arc for arc in self.path[j - 1] if arc not in self.path[j]], self.n, self.q)
                        .is_cascade()):
                    print(f"Cascade of waves fails at swell {j}")
                    return False
                else:
                    print(f"Set partition fails at swell {j}")
                    return False
        else:
            return True

    def make_beach_matrix(self):
        """
        Constructs an (n - 1) by floor(k) matrix representation of the semibeach.

        Returns:
            Matrix: The matrix representation of the semibeach.
        """
        b = self.semibeach
        mat_b = zeros(self.n - 1, floor(self.k))

        for i in range(floor(self.k)):
            for j in range(2 * i + 2, int(2 * self.k) + 1):
                if j % 2 == 1:
                    for arc in b[i][j]:
                        if arc not in b[i][j - 1] and Arc(arc).dim > -1:
                            mat_b[Arc(arc).re - 1, i] = Arc(arc).lab
                else:
                    for arc in b[i][j]:
                        if arc not in b[i][j - 1] and Arc(arc).dim > -1:
                            mat_b[Arc(arc).le - 1, i] = Arc(arc).lab

        return mat_b

    def make_beach_matrix_without_fsp(self):
        """
        Constructs an (n - 1) by floor(k) matrix representation of the semibeach without the data of the final set
        partition.

        Returns:
            numpy.ndarray: The matrix representation of the semibeach without its final set partition.
        """
        b = self.semibeach
        mat_b = zeros(self.n - 1, floor(self.k))

        for i in range(floor(self.k)):
            for j in range(2 * i + 2, int(2 * self.k) + 1):
                if j % 2 == 1:
                    for arc in b[i][j]:
                        if arc not in b[i][j - 1] and Arc(arc).dim > -1 and arc not in b[i][-1]:
                            mat_b[Arc(arc).re - 1, i] = Arc(arc).lab
                else:
                    for arc in b[i][j]:
                        if arc not in b[i][j - 1] and Arc(arc).dim > -1 and arc not in b[i][-1]:
                            mat_b[Arc(arc).le - 1, i] = Arc(arc).lab

        return mat_b
