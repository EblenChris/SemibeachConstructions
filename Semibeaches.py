from SemibeachGenerator import *
from pathlib import Path
import re


class Semibeach:
    """
    Represents a semibeach structure, which is a combinatorial object related to arc diagrams and partitions.
    """
    def __init__(self, semibeach: Tuple[List[List[List[Tuple[Tuple[int, int], Any]]]], List[Tuple]], q: int):
        """
        Initializes a Semibeach object.

        Args:
            semibeach (Tuple): A tuple containing a matrix representation of the semibeach and associated shifts.
            q (int): The order of the finite field.
        """
        self.semibeach = semibeach[0]
        self.shifts = semibeach[1]
        self.row_size = len(self.semibeach)
        self.col_size = len(self.semibeach[0])
        self.n = Arc(Wave(self.semibeach[0][2]).face).re
        self.k = float((self.col_size - 1) / 2)
        self.q = q
        self.path = [[Wave(self.semibeach[i][j]).face for i in range((j - 2) // 2 + 1) if len(self.semibeach[i][j]) > 0
                      and Arc(Wave(self.semibeach[i][j]).face).dim > -1] for j in range(self.col_size)]

    @property
    def semibeach_labels(self) -> List[List[List[Any]]]:
        """
        Computes labels for the semibeach.

        Returns:
            List: An array whose entries contains the labels associated with the corresponding entry of the semibeach.
            The list of labels in entry (i, j) corresponds to the labels of the arcs in entry (i, j) of the semibeach
            in descending arc order.
        """
        labels = [[[] for _ in range(self.col_size)] for _ in range(self.row_size)]
        for i in range(self.row_size):
            for j in range(self.col_size):
                if j == 2 * i + 3 and ((self.n, self.n), 0) in self.semibeach[i][j]:
                    labels[i][j].append(0)
                else:
                    for h in range(len(self.semibeach[i][j])):
                        labels[i][j].append(self.semibeach[i][j][h][1])
                labels[i][j].reverse()
        return labels

    def latex_preamble(self) -> str:
        """
        Generates the LaTeX preamble for the document, setting up the document class, packages, geometry,
        and other necessary configurations.

        Returns:
            str: The LaTeX preamble as a formatted string.
        """
        preamble = (f"\\documentclass[12pt, a4paper, landscape]{{article}}\n"
                    f"\\renewcommand{{\\baselinestretch}}{{1.25}}\n"
                    f"\\usepackage{{amsmath, amsthm, amssymb, array, bm, tikz, tikz-cd, bbm, mathtools, bbold, "
                    f"verbatim, booktabs, colortbl, dsfont}}\n"
                    f"\\usepackage[margin=2cm]{{geometry}}\n"
                    f"\\usetikzlibrary{{positioning}}\n"
                    f"\\setcounter{{MaxMatrixCols}}{{{max(11, self.col_size - 2)}}}\n"
                    f"\\renewcommand{{\\vec}}[1]{{\\mathbf{{#1}}}}\n"
                    f"\\begin{{document}}")
        return preamble

    def latex_semibeach(self, semibeach_name: str) -> str:
        """
        Generates the LaTeX code to visualize the semibeach structure as a table with TikZ graphics.

        Args:
            semibeach_name (str): The name of the semibeach.

        Returns:
            str: The LaTeX code representing the semibeach structure.
        """
        scale = 15 / (self.n * self.col_size)
        table_header = (f"\\begin{{center}}\n${semibeach_name} = "
                        f"$\n \\begin{{tabular}}{{|{'c' * (self.col_size - 2)}|}}\n")

        def tikz_cell(arcs):
            dots = " ".join(f"({dot}, 0) circle (2pt);" for dot in range(self.n))
            arcs_tex = "\n".join(
                f"\\draw ({arc.le - 1}, 0) to[out=60,in=120{',distance=8mm' if arc.dim == -1 else ''}] "
                f"({arc.re - 1}, 0);"
                for arc in map(Arc, arcs)
            )
            return f"$\\begin{{tikzpicture}}[scale={scale}]\\fill[black]\n{dots}\n{arcs_tex}\n\\end{{tikzpicture}}$"

        table_body = " \\\\\n".join(
            " & ".join(tikz_cell(self.semibeach[i][j]) for j in range(2, self.col_size))
            for i in range(self.row_size)
        )

        return f"{table_header}{table_body}\n\\end{{tabular}}\n\\end{{center}}\n"

    def latex_labels(self, semibeach_name: str) -> str:
        """
        Generates the LaTeX code for labeling arcs in the semibeach structure.

        Args:
            semibeach_name (str): The name of the semibeach.

        Returns:
            str: The LaTeX representation of the labels.
        """

        def format_label(label):
            """Formats individual labels for LaTeX output."""
            if label == 0 or self.q == 2:
                return str(label)
            else:
                label_str = str(label)
                return f"{label_str[0]}_{{{label_str[1:]}}}" if len(label_str) > 1 else label_str

        def format_row(row):
            """Formats an entire row of labels for LaTeX."""
            formatted_entries = []
            for j, label_list in enumerate(row):
                if j < 2:  # Skip first two columns
                    continue
                if not label_list:
                    formatted_entries.append("()")
                else:
                    formatted_entries.append(f"({', '.join(format_label(label) for label in label_list)})")
            return " & ".join(formatted_entries)

        labels_tex = f"\\[\n\\text{{Labels}}({semibeach_name}) = \n \\begin{{pmatrix}}\n"
        labels_tex += " \\\\\n".join(format_row(row) for row in self.semibeach_labels)
        labels_tex += "\n\\end{pmatrix}\n\\]\n"

        return labels_tex

    def latex_matrix(self, m: List) -> str:
        """
        Converts a given matrix into its LaTeX representation.

        Args:
            m (List): A 2D list representing the matrix.

        Returns:
            str: A LaTeX formatted string representing the matrix.
        """
        matrix_tex = "\n \\begin{pmatrix}\n"
        for i, row in enumerate(m):
            for j, entry in enumerate(row):
                if j < len(row) - 1:
                    matrix_tex += f"{str(m[i][j])[0]}_{{{str(m[i][j])[1:]}}} & " if len(str(m[i][j])[1:]) > 0 else f"{str(m[i][j])[0]} & "
                elif i < self.n - 1:
                    matrix_tex += f"{str(m[i][j])[0]}_{{{str(m[i][j])[1:]}}} \\\\ \n " if len(str(m[i][j])[1:]) > 0 else f"{str(m[i][j])[0]} \\\\ \n "
                else:
                    matrix_tex += f"{str(m[i][j])[0]}_{{{str(m[i][j])[1:]}}} \n " if len(str(m[i][j])[1:]) > 0 else f"{str(m[i][j])[0]} \n"
        matrix_tex += "\\end{pmatrix}\n\\]\n"
        return matrix_tex

    def latex_basis_elements(self, semibeach_name: str) -> str:
        """
        Generates the LaTeX representation of basis elements for a given semibeach submodule.

        Args:
            semibeach_name (str): The name of the semibeach submodule.

        Returns:
            str: A LaTeX formatted string representing the basis elements.
        """
        theta_arg, tensors = SemibeachSubmodule(self.n, self.k, self.q).basis_element([self.semibeach, self.shifts])

        def format_expression(expr):
            """Applies common formatting rules to an expression."""
            expr = re.sub(r"f(\d+)", r"f_{\1}", expr)
            expr = re.sub(r"([st])(\d+)\[(\d+), (\d+)]",
                          lambda m: f"({m.group(1)}_{{{int(m.group(2)) + 1}}})_{{{int(m.group(3)) + 1}}}", expr)
            expr = expr.replace("*", "")
            expr = re.sub(r"(\d+)\.0+", r"\1", expr)
            expr = re.sub(r"/([^/()]+)", r"\1^{-1}", expr)
            expr = re.sub(r"/\(([^()]+)\)", r"(\1)^{-1}", expr)
            expr = re.sub(r"\b1\*", "", expr)
            expr = re.sub(r"\*\*", "^", expr)
            expr = re.sub(r"\^(-?\d+)\^{-1}|\^{-1}\^(-?\d+)",
                          lambda m: f"^{{-{m.group(1) or m.group(2)}}}", expr)
            return expr

        formatted_theta_arg = format_expression(str(theta_arg))

        def format_tensor(matrix):
            """Formats a tensor matrix as a LaTeX pmatrix."""
            elements = [format_expression(str(row[-1])) for row in matrix.tolist()]
            return "\n\\begin{pmatrix}\n" + " \\\\\n".join(elements) + "\n\\end{pmatrix}"

        formatted_tensors = " \\otimes ".join(map(format_tensor, tensors)) + " \\otimes \\mathds 1"

        return (f"\\begin{{align*}}\n\\nu^{{{semibeach_name}}}(\\vec t) = "
                f"\\sum_{{\\vec s}} \\vartheta({formatted_theta_arg}) {formatted_tensors}\n\\end{{align*}}")

    def tex_draw(self, file_name: str, second_beach: Any) -> None:
        """
        Generates a LaTeX file visualizing the semibeach structure, its associated module, and the associated beach
        matrix when the semibeach is a beach.

        Args:
            file_name (str): The name of the file to be created.
            second_beach (Any): The second beach to construct in LaTeX.
        """
        # creates a folder called file_name to store the file file_name.txt
        output_file = Path(f"{file_name}/{file_name}.tex")
        output_file.parent.mkdir(exist_ok=True, parents=True)
        with open(output_file, "w") as file:
            file.write(self.latex_preamble())
            if second_beach is None:
                file.write(self.latex_semibeach("B"))
                file.write(self.latex_labels("B"))
                file.write("\\[\n\\mathrm{Mat}(B) = ")
                file.write(self.latex_matrix(self.make_beach_matrix().tolist()))
                file.write(self.latex_basis_elements("B"))
            else:
                file.write(self.latex_semibeach("B_1"))
                file.write(self.latex_labels("B_1"))
                file.write(Semibeach(second_beach, self.q).latex_semibeach("B_2"))
                file.write(Semibeach(second_beach, self.q).latex_labels("B_2"))
                file.write("\\newpage")
                file.write("\\[\n(\\overline{\\mathrm{Mat}(B_1)} \\mid \\mathrm{Mat}(B_2)) = ")
                mat_1 = self.make_beach_matrix_without_fsp().tolist()
                mat_2 = Semibeach(second_beach, self.q).make_beach_matrix().tolist()
                mat = [row1 + row2 for row1, row2 in zip(mat_1, mat_2)]
                file.write(self.latex_matrix(mat))
            file.write("\\end{document}\n")

    def cascade_and_partition_checker(self) -> bool:
        """
        Checks whether the structure satisfies the cascade and set partition conditions. Prints diagnostic messages
        indicating where the conditions fail.

        Returns:
            bool: True if every two successive swells" symmetric difference is a cascade of waves and every swell"s
            union of nontrivial faces is a set partition, False otherwise.
        """
        for j in range(1, self.col_size):
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
        mat_b = zeros(self.n - 1, self.row_size)

        for i in range(self.row_size):
            for j in range(2 * i + 2, self.col_size):
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
        mat_b = zeros(self.n - 1, self.row_size)

        for i in range(self.row_size):
            for j in range(2 * i + 2, self.col_size):
                if j % 2 == 1:
                    for arc in b[i][j]:
                        if arc not in b[i][j - 1] and Arc(arc).dim > -1 and arc not in b[i][-1]:
                            mat_b[Arc(arc).re - 1, i] = Arc(arc).lab
                else:
                    for arc in b[i][j]:
                        if arc not in b[i][j - 1] and Arc(arc).dim > -1 and arc not in b[i][-1]:
                            mat_b[Arc(arc).le - 1, i] = Arc(arc).lab

        return mat_b


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

    def basis_element(self, semibeach: Tuple[List[List[Any]], List[Tuple]]) -> Tuple[Any, List[Any]]:
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
                curr_t_part (List[Matrix]): The current \mathcal T^B matrices.
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
                l = min(n - 1, Arc(seed_arc_left).le) - 1
                curr_tensors += [curr_t_part[0] * proj(s_mat(0), 0, l)]
                curr_arg += - Arc(seed_arc_left).lab * s_mat(0)[l, n - 1]
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
                l = min(n - 1, Arc(seed_arc_left).le) - 1
                curr_tensors += [curr_t_part[num_rows] * proj(s_mat(num_rows), 0, l)]
                curr_arg += - Arc(seed_arc_left).lab * s_mat(num_rows)[l, n - 1]
                new_t_part = curr_t_part[0:num_rows]
                return recursion_basis_elements(curr_arg, curr_tensors, new_t_part, m - 1)
            elif len(c) == 3 and curr_swell % 2 == 0:
                losing_arc, resolving_arc = Shift(c).losing_slice()
                winning_arc, winning_row = Shift(c).winning_arc(), Shift(c).winning_row()
                l = min(Arc(resolving_arc).re - 1, Arc(resolving_arc).le) - 1
                curr_arg += - Arc(resolving_arc).lab * s_mat(losing_row)[l, n - 1]
                t_los = simplify((proj(s_mat(losing_row), Arc(losing_arc).le, l) * curr_t_part[losing_row]))
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
                l = Arc(resolving_arc).le - 1
                curr_arg += (Arc(losing_arc).lab *
                             curr_t_part[winning_row][l, n - 1])
                t_los = simplify(scale(curr_t_part[Shift(c).losing_row()], (1 / Arc(losing_arc).lab) *
                                       Arc(resolving_arc).lab) * proj(curr_t_part[winning_row], l + 1, n - 1))
                new_t_part = []
                for h in range(num_rows + 1):
                    if h != losing_row:
                        new_t_part.append(curr_t_part[h])
                    else:
                        new_t_part.append(t_los)
                return recursion_basis_elements(curr_arg, curr_tensors, new_t_part, m - 1)

        return recursion_basis_elements(0, [], t_b, len(shifts) - 1)
