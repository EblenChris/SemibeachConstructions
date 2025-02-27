from ArcDiagrams import *
from math import floor


class SemibeachSubmodule:
    """
    Represents the semibeach submodule.
    """
    def __init__(self, n: int, k: float, q: int):
        """
        Initializes a SemibeachSubmodule instance.

        Args:
            n (int): The parameter n for UT_n(q).
            k (float): The parameter k for the level of the Bratteli diagram.
            q (int): The order of the finite field.
        """
        self.n = n
        self.k = k
        self.q = q
        self.row_size = floor(k)
        self.col_size = int(2 * k + 1)
        self.s_symbols = list(symbols("s0:%d" % self.row_size))
        self.t_symbols = list(symbols("t0:%d" % self.row_size))

    def basis_element(self, semibeach: Tuple[List[List[List[Tuple[Tuple[int, int], Any]]]], List[Tuple]]) -> (
            Tuple)[Any, List[Any]]:
        """
        Computes the structure of the basis elements of the submodule.

        Args:
            semibeach: A tuple consisting of a semibeach and its shifts.

        Returns:
             A tuple containing the argument for the nontrivial group homomorphism theta and the form of the tensor
             factors for basis elements in the given semibeach module.
        """
        n = self.n
        b, shifts = semibeach

        def crash_right(row: int) -> bool:
            for shift in shifts:
                if all([Arc(Wave(Shift(shift).losing_slice()).face).dim == -1, Shift(shift).swell() % 2 == 1,
                        Shift(shift).losing_row() == row]):
                    return True
            return False

        def proj(v: Matrix, m: int, r: int) -> Matrix:
            v_0 = eye(self.n)
            for h in range(m, r + 1):
                v_0[h, self.n - 1] = v[h, self.n - 1]
            return v_0

        def scale(v: Matrix, a: Any) -> Matrix:
            v_0 = eye(self.n)
            for h in range(self.n - 1):
                v_0[h, n - 1] = a * v[h, n - 1]
            return v_0

        def t_mat(h: int) -> Matrix:
            return proj(Matrix(MatrixSymbol(self.t_symbols[h], n, n)), 0, n - 2)

        def s_mat(h: int) -> Matrix:
            return proj(Matrix(MatrixSymbol(self.s_symbols[h], n, n)), 0, n - 2)

        t_b = []

        for i in range(self.row_size):
            if crash_right(i):
                t_b.append(eye(n))
            else:
                right = Arc(Wave(b[i][self.col_size - 1]).face).re
                left = min(Arc(Wave(b[i][self.col_size - 1]).face).le, right - 1)
                t = proj(t_mat(i), left, right - 1)
                t[right - 1, n - 1] = 1
                t_b.append(t)

        def recursion_basis_elements(curr_arg: Any, curr_tensors: List[Matrix], curr_t_part: List[Matrix], m: int):
            """
            Recursive function to construct basis elements.

            Args:
                curr_arg (Any): The current argument of the nontrivial homomorphism theta.
                curr_tensors (List[Matrix]): The list of current tensor factors of basis elements.
                curr_t_part (List[Matrix]): The current \\mathcal T^B matrices.
                m (int): The current index in the recursion.
            Returns:
                Updated argument and tensors.
            """
            c = shifts[m]
            curr_swell = Shift(c).swell()
            num_rows = (curr_swell - 2) // 2
            losing_row = Shift(c).losing_row()
            if m == 0:
                # define submodule for initial left seed and return basis elements
                seed_arc_left = Wave(Shift(shifts[0]).losing_slice()).face
                le = min(n - 1, Arc(seed_arc_left).le) - 1
                curr_tensors += [curr_t_part[0] * proj(s_mat(0), 0, le)]
                curr_arg += - Arc(seed_arc_left).lab * s_mat(0)[le, n - 1]
                return curr_arg, curr_tensors
            elif len(c) == 2 and curr_swell % 2 == 1:
                seed_arc_left, seed_arc_right = Shift(c).losing_slice()
                seed_arc_left_label = Arc(seed_arc_left).lab if (Arc(seed_arc_left).lab != 0) else 1
                new_t_part = (curr_t_part[0:num_rows] + [simplify(scale(curr_t_part[num_rows],
                                                                        (1 / seed_arc_left_label) *
                                                                        Arc(seed_arc_right).lab))])
                return recursion_basis_elements(curr_arg, curr_tensors, new_t_part, m - 1)
            elif len(c) == 2 and curr_swell % 2 == 0:
                seed_arc_left = Wave(Shift(c).losing_slice()).face
                le = min(n - 1, Arc(seed_arc_left).le) - 1
                curr_tensors += [curr_t_part[num_rows] * proj(s_mat(num_rows), 0, le)]
                curr_arg += - Arc(seed_arc_left).lab * s_mat(num_rows)[le, n - 1]
                new_t_part = curr_t_part[0:num_rows]
                return recursion_basis_elements(curr_arg, curr_tensors, new_t_part, m - 1)
            elif len(c) == 3 and curr_swell % 2 == 0:
                losing_arc, resolving_arc = Shift(c).losing_slice()
                winning_arc, winning_row = Shift(c).winning_arc(), Shift(c).winning_row()
                le = min(Arc(resolving_arc).re - 1, Arc(resolving_arc).le) - 1
                curr_arg += - Arc(resolving_arc).lab * s_mat(losing_row)[le, n - 1]
                t_los = simplify((proj(s_mat(losing_row), Arc(losing_arc).le, le) * curr_t_part[losing_row]))
                t_win = simplify(curr_t_part[winning_row] * scale(t_los, - (1 / Arc(winning_arc).lab)
                                                                  * Arc(losing_arc).lab))
                new_t_part = []
                for h in range(num_rows + 1):
                    if h not in {losing_row, winning_row}:
                        new_t_part.append(curr_t_part[h])
                    elif h == losing_row:
                        new_t_part.append(t_los)
                    else:
                        new_t_part.append(t_win)
                return recursion_basis_elements(curr_arg, curr_tensors, new_t_part, m - 1)
            elif len(c) == 3 and curr_swell % 2 == 1:
                losing_arc, resolving_arc = Shift(c).losing_slice()
                winning_row = Shift(c).winning_row()
                le = Arc(resolving_arc).le - 1
                curr_arg += (Arc(losing_arc).lab *
                             curr_t_part[winning_row][le, n - 1])
                t_los = simplify(scale(curr_t_part[Shift(c).losing_row()], (1 / Arc(losing_arc).lab) *
                                       Arc(resolving_arc).lab) * proj(curr_t_part[winning_row], le + 1, n - 1))
                new_t_part = []
                for h in range(num_rows + 1):
                    if h != losing_row:
                        new_t_part.append(curr_t_part[h])
                    else:
                        new_t_part.append(t_los)
                return recursion_basis_elements(curr_arg, curr_tensors, new_t_part, m - 1)

        return recursion_basis_elements(0, [], t_b, len(shifts) - 1)


class Shift:
    """
    Represents a shift of the semibeach.
    """
    def __init__(self, shift: Tuple):
        """
        Initializes a Shift instance.

        Args:
            shift: A tuple representing the shift data.
        """
        self.shift = shift

    def swell(self) -> int:
        """
        Returns the swell value of the shift.

        Returns:
             Swell value as an integer.
        """
        return self.shift[2] if len(self.shift) == 3 else self.shift[1]

    def winning_arc(self) -> Tuple[Tuple[int, int], Any]:
        """
        Returns the winning arc of the shift.

        Returns:
             Winning arc.
        """
        return self.shift[0][0]

    def winning_row(self) -> int:
        """
        Returns the winning row of the shift.

        Returns:
            Winning row index (0-based).
        """
        return int(self.shift[0][1])

    def losing_slice(self) -> Tuple[Tuple[int, int], Any]:
        """
        Returns the losing slice of the shift.

        Returns:
             Losing slice.
        """
        if len(self.shift) == 2 and Shift(self.shift).swell() % 2 == 0:
            return [self.shift[0]]
        elif len(self.shift) == 2 and Shift(self.shift).swell() % 2 == 1:
            return self.shift[0]
        if len(self.shift) == 3:
            return self.shift[1][0]

    def losing_row(self) -> int:
        """
        Returns the losing row of the shift.

        Returns:
            Losing row index (0-based).
        """
        return (self.shift[1] // 2) - 1 if len(self.shift) == 2 else self.shift[1][1]
