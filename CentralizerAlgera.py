import itertools
from LaTeX import *
from datetime import datetime
import pandas as pd


class CentralizerAlgebra:
    """
    A class representing the centralizer algebra for UT_n(q) or UT_{n - 1}(q) at level k. This class provides methods
    for computing the dimension of the algebra, counting beach maps, and generating related data such as cascades,
    paths, and coefficients. Additionally, it facilitates exporting results to spreadsheets and analyzing formula
    coefficients.

    Attributes:
        n (int): The size of the unitriangular group.
        k (float): The parameter k for the level of the Bratteli diagram.
        q (int): The parameter q used for UT_n(q).
        dimension (int): The dimension of the centralizer algebra.

    Properties:
        set_partitions_of_n (tuple): All (unlabeled) set partitions of {1, 2, ..., n}.
        set_partitions_of_n_minus_1 (tuple): All (unlabeled) set partitions of {1, 2, ..., n - 1}.
        cascades (list): All valid cascades of waves for the given Args.

    Methods:
        num_beach_maps_and_beaches(extra_q_values: List[int]) -> Tuple:
            Computes the number of beach maps and the number of beaches for a list of additional q values.

        get_paths() -> list:
            Generates all valid paths of length k in the Bratteli diagram for set partitions of n.

        get_spreadsheet(extra_n_values: List[int], extra_k_values: List[float], extra_q_values: List[int]):
            Generates a spreadsheet of results including dimensions, number of beach maps, and ratios.

        get_difference_coeffs() -> None:
            Computes the coefficients for binomial terms in a formula describing the difference between the dimension
            of the algebra and the number of beach maps.
    """
    def __init__(self, n: int, k: float, q: int):
        """
        Initializes the CentralizerAlgebra object with the given Args n, k, and q.

        Args:
            n (int): The parameter n for UT_n(q).
            k (float): The parameter k for the length of the Bratteli diagram.
            q (int): The parameter q used for labeling arcs.
        """
        self.n = n
        self.k = k
        self.q = q
        self.dimension = sum(comb(n - 1, i) * prod([q ** (int(2 * k) - j) - 1 for j in range(1, i + 1)])
                             for i in range(min(int(2 * k) - 1, n - 1) + 1))

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

    @property
    @lru_cache(maxsize=None)
    def cascades(self) -> List[set]:
        """
        Recursively generates all cascades of waves for the given Args.

        Returns:
            list: A list of sets, where each set represents a cascade of waves.
        """
        def generate_cascades(cascades, forming_cascades):
            still_forming_cascades = []

            for forming_cascade in forming_cascades:
                left, right = Arc(forming_cascade[-1]).le, Arc(forming_cascade[-1]).re

                if (len(forming_cascade) - 1) % 2 == 0:
                    for arc_choice in [((left, j), 1) for j in range(left + 1, right)]:
                        still_forming_cascades.append(forming_cascade.copy() + [arc_choice])
                else:
                    for arc_choice in [((j, right), 1) for j in range(left + 1, right)]:
                        still_forming_cascades.append(forming_cascade.copy() + [arc_choice])

            if not still_forming_cascades:
                return cascades + forming_cascades

            return generate_cascades(cascades + forming_cascades, still_forming_cascades)

        initial_cascades = [[((i, self.n), 1)] for i in range(1, self.n)]
        return [set()] + [set(cascade) for cascade in generate_cascades([], initial_cascades)]

    def num_beach_maps_and_beaches(self, extra_q_values: List[int]) -> Tuple[List[int], List[int]]:
        """
        Computes the number of beach maps and the number of beaches for a list of additional q values.

        Args:
            extra_q_values (List[int]): A list of additional q values.

        Returns:
            Tuple[list[int], list[int]]: A tuple containing two lists:
                - The number of beach maps for each q value.
                - The number of beaches for each q value.
        """
        q_values = [self.q] + extra_q_values
        cascades = self.cascades

        starting_weight_dict = {tuple([((i, self.n), 1)]): [[1 for _ in range(len(q_values))]]
                                for i in range(1, self.n)}
        starting_weight_dict[tuple()] = [[1 for _ in range(len(q_values))]]

        def helper(curr_weight_dict, depth):
            if depth == self.k:
                print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}: n = {self.n}, '
                      f'computing number of beaches and beach maps')

                num_beach_maps, num_beaches = [0 for _ in range(len(q_values))], [0 for _ in range(len(q_values))]

                for fsp in curr_weight_dict:
                    num_beaches_to_fsp = [sum([curr_weight_dict[tuple(fsp)][j][i]
                                               for j in range(len(curr_weight_dict[tuple(fsp)]))])
                                          for i in range(len(q_values))]
                    for i, q in enumerate(q_values):
                        num_beach_maps[i] += ((q - 1) ** len(fsp)) * (num_beaches_to_fsp[i] ** 2)
                        num_beaches[i] += ((q - 1) ** len(fsp)) * num_beaches_to_fsp[i]

                return num_beach_maps, num_beaches

            else:
                print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}: n = {self.n}, recursion step k = {depth + 0.5}')
                if int(2 * depth) % 2 == 1:
                    next_fsps = [sp for sp in self.set_partitions_of_n if len(sp) <= int(depth + 0.5) and
                                 (len(sp) != int(depth + 0.5) or self.n in ArcDiagram(sp, self.n, 2).rights)]
                else:
                    next_fsps = [sp for sp in self.set_partitions_of_n_minus_1 if len(sp) <= floor(depth + 0.5)]

                next_weight_dict = {tuple(next_fsp): [] for next_fsp in next_fsps}
                for j, sp in enumerate(curr_weight_dict.keys()):
                    if (j + 1) % 1000 == 0:
                        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}: n = {self.n}, '
                              f'recursion step k = {depth + 0.5}, getting weights for new edges from {j + 1}th '
                              f'set partition of {len(curr_weight_dict.keys())}')
                    for next_fsp in next_fsps:
                        curr_minus_next = set(sp) - set(next_fsp)
                        next_minus_curr = set(next_fsp) - set(sp)
                        if curr_minus_next.union(next_minus_curr) in cascades:
                            curr_and_next = set(sp) & set(next_fsp)
                            if int(2 * depth) % 2 == 1:
                                curr_minus_next_ad = ArcDiagram(list(curr_minus_next), self.n, 2)
                                next_minus_curr_ad = ArcDiagram(list(next_minus_curr), self.n, 2)
                                next_weight = [((q ** (len(curr_minus_next_ad.crossings(list(curr_and_next))) - len(
                                    next_minus_curr_ad.crossings(list(curr_and_next))))) * (
                                             (q - 1) ** len(list(curr_minus_next)))) for i, q in enumerate(q_values)]
                                for curr_path_weight in curr_weight_dict[sp]:
                                    next_weight_dict[tuple(next_fsp)].append([curr_path_weight[i] * next_weight[i]
                                                                              for i in range(len(q_values))])
                            else:
                                curr_and_next_ad = ArcDiagram(list(curr_and_next), self.n, 2)
                                next_weight = [((q ** (len(curr_and_next_ad.crossings(list(curr_minus_next))) - len(
                                    curr_and_next_ad.crossings(list(next_minus_curr))))) * (
                                        (q - 1) ** len(list(curr_minus_next)))) for i, q in enumerate(q_values)]
                                for curr_path_weight in curr_weight_dict[sp]:
                                    next_weight_dict[tuple(next_fsp)].append([curr_path_weight[i] * next_weight[i]
                                                                              for i in range(len(q_values))])

                return helper(next_weight_dict, depth + 0.5)

        return helper(starting_weight_dict, 1)

    @lru_cache(maxsize=None)
    def get_paths(self) -> List[List[List[Tuple[Tuple[int, int], Any]]]]:
        """
        Generates all valid paths in the Bratteli diagram up to the given level k.

        Paths are sequences of set partitions connected by edges.

        Returns:
            List[List[List[Tuple[Tuple[int, int], Any]]]]: A list of paths, where each path is represented as a list
            of set partitions.
        """
        cascades = self.cascades

        def get_paths_helper(curr_paths, depth):
            if depth == self.k:
                return curr_paths
            else:
                print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}: n = {self.n}, recursion step k = {depth + 0.5}')
                if int(2 * (depth + 0.5)) % 2 == 0:
                    next_superchars = [set_partition for set_partition in self.set_partitions_of_n
                                       if len(set_partition) <= floor(depth + 0.5)]
                else:
                    next_superchars = [set_partition for set_partition in self.set_partitions_of_n_minus_1
                                       if len(set_partition) <= floor(depth + 0.5)]
                next_paths = []
                for path in curr_paths:
                    for superchar in next_superchars:
                        if (set(path[-1]) - set(superchar)).union(set(superchar) - set(path[-1])) in cascades:
                            next_paths.append(path + [superchar])
                return get_paths_helper(next_paths, depth + 0.5)

        return get_paths_helper([[[]]] + [[[((i, self.n), 1)]] for i in range(1, self.n)], 1)

    def get_spreadsheet(self, extra_n: List[int], extra_k: List[float], extra_q: List[int]) -> None:
        """
        Generates a spreadsheet containing results such as dimensions, the number of beach maps,
        the number of beaches, and their ratios for various parameter values.

        Args:
            extra_n (List[int]): Additional n values to include in the computation.
            extra_k (List[float]): Additional k values to include in the computation.
            extra_q (List[int]): Additional q values to include in the computation.

        Returns:
            None: Saves the result to a CSV file.
        """
        n_values, k_values, q_values = [self.n] + extra_n, [self.k] + extra_k, [self.q] + extra_q
        res = []

        for n, k in itertools.product(n_values, k_values):
            num_beach_maps, num_beaches = (CentralizerAlgebra(n, k, self.q).num_beach_maps_and_beaches(extra_q))
            for i, q in enumerate(q_values):
                dimension = CentralizerAlgebra(n, k, q).dimension
                beach_maps, beaches = num_beach_maps[i], num_beaches[i]

                res.append([n, k, q, dimension, beach_maps, beaches, beach_maps / dimension])

        df = pd.DataFrame(res, columns=["n", "k", "q", "dimension", "beach maps", "beaches", "ratio"])

        df.to_csv("centralizer_algebra_results.csv", index=False)

    def is_rothe_diagram(self, r: Matrix) -> bool:
        """
        Checks whether a given matrix represents a Rothe diagram.

        A Rothe diagram is a combinatorial object that encodes the positions of inversions in a permutation.
        This function verifies the Rothe diagram condition by ensuring that each nonzero entry satisfies
        specific row and column constraints.

        Args:
            r (Matrix): A matrix representation of the diagram.

        Returns:
            bool: True if the matrix satisfies the Rothe diagram conditions, False otherwise.
                  If the condition fails, it prints the failing entry's position.
        """
        for j in range(int(2 * self.k - 1)):
            for i in range(self.n - 1):
                if r[self.n - i - 2, j] != 0 and {r[self.n - i - 2, h] for h in range(j)}.union({0}) == {0}:
                    if {r[h, j] for h in range(self.n - i - 2)}.union({0}) != {0}:
                        print(f'Fails at entry {self.n - i - 2}, {j}')
                        return False
        return True

    def get_sheets_formula(self, file_name: str, coeffs: List[List[int]]) -> None:
        offset_a = 3 if int(2 * self.k) % 2 == 0 else 4
        offset_c = 5 if int(2 * self.k) % 2 == 0 else 7

        sheets_formula = '= ' + ' + '.join(
            '(' + ' + '.join(
                f'{coeff} * COMBIN(A2 - 1, {offset_a + j})' for j, coeff in enumerate(row)
            ) + f') * (C2 - 1)^{offset_c + i}'
            for i, row in enumerate(coeffs)
        )

        output_file = Path(f"{file_name}.txt")
        with open(output_file, "w") as file:
            file.write(sheets_formula)

    def get_difference_coeffs(self) -> Tuple[List[List[int]], str, str]:
        """
       Computes the coefficients of binomial terms for a formula describing the difference between
       the dimension of the algebra and the number of beach maps.

       The entry self.get_difference_coeffs()[0][i][j] is the coefficient of \binom{n - 1}{j}(q - 1)^i in the formula
       for the difference between the dimension of the centralizer algebra and the number of beach maps for n
       and k.

       Returns:
           Tuple(List[List[int]], str, str): A tuple containing:
            - A list of the coefficients A_{i, j}(n, k)
            - The difference formula in Google sheets format for easy data verification.
            - The difference formula in LaTeX format.
       """

        extra_n = [n for n in range(int(2 * self.k) + 1, int(4 * self.k - 3))] if int(2 * self.k) % 2 == 0 \
            else [n for n in range(int(2 * self.k) + 1, int(4 * self.k - 4))]

        extra_q = [q for q in range(self.q + 1, comb(int(2 * self.k), 2) - 5 + self.q)] \
            if int(2 * self.k) % 2 == 0 \
            else [q for q in range(self.q + 1, comb(int(2 * self.k), 2) - 6 + self.q)]

        dimensions = [[CentralizerAlgebra(n, self.k, q).dimension for q in [self.q] + extra_q]
                      for n in [self.n] + extra_n]

        num_beach_maps, num_beaches = ([0 for _ in range(len([self.n] + extra_n))],
                                       [0 for _ in range(len([self.n] + extra_n))])

        for i, n in enumerate([self.n] + extra_n):
            print(f'Computing n = {n} of {([self.n] + extra_n)[-1]}')
            num_beach_maps[i], num_beaches[i] = (CentralizerAlgebra(n, self.k, self.q)
                                                 .num_beach_maps_and_beaches(extra_q))

        # matrix of (q - 1)^i values for the various q
        q_mat = Matrix([[(q - 1) ** i for i in range(5, comb(int(2 * self.k), 2))]
                        for q in [self.q] + extra_q]) if int(2 * self.k) % 2 == 0 \
            else Matrix([[(q - 1) ** i for i in range(7, comb(int(2 * self.k), 2))]
                        for q in [self.q] + extra_q])

        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}: Computing coefficients')

        difference_mats = [Matrix([[dimensions[j][i] - num_beach_maps[j][i]]
                                   for i in range(len([self.q] + extra_q))])
                           for j in range(len([self.n] + extra_n))]

        sols = [q_mat.LUsolve(difference_mats[j]) for j in range(len([self.n] + extra_n))]

        binomial_mat = Matrix([[comb(n - 1, i) for i in range(3, int(2 * self.k))]
                               for n in [self.n] + extra_n]) if int(2 * self.k) % 2 == 0 \
            else Matrix([[comb(n - 1, i) for i in range(4, int(2 * self.k))] for n in [self.n] + extra_n])

        coeffs = []
        if int(2 * self.k) % 2 == 0:
            for i in range(comb(int(2 * self.k), 2) - 5):
                coeffs.append(binomial_mat.LUsolve(Matrix([sols[j][i] for j in range(len([self.n] + extra_n))])))
        else:
            for i in range(comb(int(2 * self.k), 2) - 7):
                coeffs.append(binomial_mat.LUsolve(Matrix([sols[j][i] for j in range(len([self.n] + extra_n))])))

        # generate a spreadsheet with the results for q = 2 for reference and verification of the generated coefficients
        results = []
        for i, n in enumerate([self.n] + extra_n):
            results.append(
                [n, self.k, 2, dimensions[i][0], num_beach_maps[i][0], num_beaches[i][0],
                 num_beach_maps[i][0] / dimensions[i][0]])

        # Convert the results to a pandas DataFrame
        df = pd.DataFrame(results, columns=["n", "k", "q", "dimension", "beach maps", "beaches", "ratio"])

        # Save the DataFrame to a CSV file (spreadsheet)
        df.to_csv("centralizer_algebra_results.csv", index=False)

        print(f'Finished at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
        return ([binom_coeff.tolist() for binom_coeff in coeffs],
                self.get_sheets_formula('difference_formula_sheets', coeffs),
                LaTeX(self.n, self.k, self.q).get_latex_formula('difference_coefficients_formula', coeffs),
                LaTeX(self.n, self.k, self.q).get_latex_table('difference_coefficients_table', coeffs))
