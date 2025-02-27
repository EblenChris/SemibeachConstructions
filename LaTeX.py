from SemibeachSubmodule import SemibeachSubmodule
from Semibeaches import *
import re
from pathlib import Path
from typing import *
from math import floor


class LaTeX:

    def __init__(self, n: int, k: float, q: int):
        """
        Initializes the LaTeX object with parameters for generating LaTeX code related to semibeach structures.

        Args:
            n (int): The size of the unitriangular group.
            k (float): The parameter k for the level of the Bratteli diagram.
            q (int): The size of the finite field.
        """
        self.n = n
        self.k = k
        self.q = q

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
                    f"\\setcounter{{MaxMatrixCols}}{{{max(11, int(2 * self.k))}}}\n"
                    f"\\renewcommand{{\\vec}}[1]{{\\mathbf{{#1}}}}\n"
                    f"\\begin{{document}}")
        return preamble

    def latex_semibeach(self, semibeach: Tuple[List[List[List[Tuple[Tuple[int, int], Any]]]], List[Tuple]],
                        semibeach_name: str) -> str:
        """
        Generates the LaTeX code to visualize the semibeach structure as a table with TikZ graphics.

        Args:
            semibeach (Tuple): The semibeach to render in LaTeX.
            semibeach_name (str): The name of the semibeach.

        Returns:
            str: The LaTeX code representing the semibeach structure.
        """
        scale = 15 / (self.n * int(2 * self.k) + 1)
        table_header = (f"\\begin{{center}}\n${semibeach_name} = "
                        f"$\n \\begin{{tabular}}{{|{'c' * (int(2 * self.k) - 1)}|}}\n")

        def tikz_cell(arcs):
            dots = " ".join(f"({dot}, 0) circle (2pt)" + (";" if dot == self.n - 1 else "") for dot in range(self.n))
            arcs_tex = "\n".join(
                f"\\draw ({arc.le - 1}, 0) to[out=60,in=120{',distance=8mm' if arc.dim == -1 else ''}] "
                f"({arc.re - 1}, 0);"
                for arc in map(Arc, arcs)
            )
            return f"$\\begin{{tikzpicture}}[scale={scale}]\\fill[black]\n{dots}\n{arcs_tex}\n\\end{{tikzpicture}}$"

        table_body = " \\\\\n".join(
            " & ".join(tikz_cell(Semibeach(semibeach, self.n, self.k, self.q).semibeach[i][j])
                       for j in range(2, int(2 * self.k) + 1))
            for i in range(floor(self.k))
        )

        return f"{table_header}{table_body}\n\\end{{tabular}}\n\\end{{center}}\n"

    def latex_labels(self, semibeach_labels: List[List[List[Any]]], semibeach_name: str) -> str:
        """
        Generates the LaTeX code for labeling arcs in the semibeach structure.

        Args:
            semibeach_labels (Tuple): The semibeach labels to render in LaTeX.
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
        labels_tex += " \\\\\n".join(format_row(row) for row in semibeach_labels)
        labels_tex += "\n\\end{pmatrix}\n\\]\n"

        return labels_tex

    @staticmethod
    def latex_matrix(m: List) -> str:
        """
        Converts a given matrix into its LaTeX representation.

        Args:
            m (List): A 2D list representing the matrix.

        Returns:
            str: A LaTeX formatted string representing the matrix.
        """
        matrix_tex = "\n \\begin{pmatrix}\n"
        for i, row in enumerate(m):
            matrix_tex += " & ".join(
                f"{str(entry)[0]}_{{{str(entry)[1:]}}}" if len(str(entry)) > 1 else f"{str(entry)[0]}"
                for entry in row
            )
            matrix_tex += " \\\\ \n" if i < len(m) - 1 else "\n"
        matrix_tex += "\\end{pmatrix}\n\\]\n"
        return matrix_tex

    def latex_basis_elements(self, semibeach: Tuple[List[List[List[Tuple[Tuple[int, int], Any]]]], List[Tuple]],
                             semibeach_name: str) -> str:
        """
        Generates the LaTeX representation of basis elements for a given semibeach submodule.

        Args:
            semibeach (Tuple): The semibeach to generate basis elements for in LaTeX.
            semibeach_name (str): The name of the semibeach submodule.

        Returns:
            str: A LaTeX formatted string representing the basis elements.
        """

        theta_arg, tensors = SemibeachSubmodule(self.n, self.k, self.q).basis_element(semibeach)

        def format_expression(expr):
            """Applies common formatting rules to an expression."""
            expr = re.sub(r"f(\d+)", r"f_{\1}", expr)
            expr = re.sub(r"([st])(\d+)\[(\d+), (\d+)]",
                          lambda m: f"({m.group(1)}_{{{int(m.group(2)) + 1}}})_{{{int(m.group(3)) + 1}}}", expr)
            expr = re.sub(r"(\d+)\.0+", r"\1", expr)
            expr = re.sub(r"\b1\s*\*", "", expr)
            expr = re.sub(r"\b-1\s*\*", "-", expr)
            expr = expr.replace("*", "")
            expr = re.sub(r"/([a-zA-Z_{}0-9]+)", r"(\1)^{-1}", expr)
            expr = re.sub(r"/\(([^()]+)\)", r"(\1)^{-1}", expr)
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

    def tex_draw(self, file_name: str, first_beach: Tuple[List[List[List[Tuple[Tuple[int, int], Any]]]], List[Tuple]],
                 second_beach: Any) -> None:
        """
        Generates a LaTeX file visualizing the semibeach structure, its associated module, and the associated beach
        matrix when the semibeach is a beach.

        Args:
            file_name (str): The name of the file to be created.
            first_beach (Tuple): The first beach to construct in LaTeX.
            second_beach (Any): The second beach to construct in LaTeX.
        """
        # creates a folder called file_name to store the file file_name.txt
        output_file = Path(f"{file_name}/{file_name}.tex")
        output_file.parent.mkdir(exist_ok=True, parents=True)
        with open(output_file, "w") as file:
            file.write(LaTeX(self.n, self.k, self.q).latex_preamble())
            if second_beach is None:
                file.write(LaTeX(self.n, self.k, self.q).latex_semibeach(first_beach, "B"))
                file.write(LaTeX(self.n, self.k, self.q).latex_labels(
                    Semibeach(first_beach, self.n, self.k, self.q).semibeach_labels, "B"))
                file.write("\\[\n\\mathrm{Mat}(B) = ")
                file.write(LaTeX(self.n, self.k, self.q).latex_matrix(
                    Semibeach(first_beach, self.n, self.k, self.q).make_beach_matrix().tolist()))
                file.write(LaTeX(self.n, self.k, self.q).latex_basis_elements((
                    Semibeach(first_beach, self.n, self.k, self.q).semibeach,
                    Semibeach(first_beach, self.n, self.k, self.q).shifts), "B"))
            else:
                file.write(LaTeX(self.n, self.k, self.q).latex_semibeach(first_beach, "B_1"))
                file.write(LaTeX(self.n, self.k, self.q).latex_labels(
                    Semibeach(first_beach, self.n, self.k, self.q).semibeach_labels, "B_1"))
                file.write(LaTeX(self.n, self.k, self.q).latex_semibeach(second_beach, "B_2"))
                file.write(LaTeX(self.n, self.k, self.q).latex_labels(
                    Semibeach(second_beach, self.n, self.k, self.q).semibeach_labels, "B_2"))
                file.write("\\newpage")
                file.write("\\[\n(\\overline{\\mathrm{Mat}(B_1)} \\mid \\mathrm{Mat}(B_2)) = ")
                mat_1 = Semibeach(first_beach, self.n, self.k, self.q).make_beach_matrix_without_fsp().tolist()
                mat_2 = Semibeach(second_beach, self.n, self.k, self.q).make_beach_matrix().tolist()
                mat = [row1 + row2 for row1, row2 in zip(mat_1, mat_2)]
                file.write(LaTeX(self.n, self.k, self.q).latex_matrix(mat))
            file.write("\\end{document}\n")

    def get_latex_formula(self, file_name: str, coeffs: List[List[int]]) -> None:
        """
        Generates a LaTeX file containing the formula for the difference of the dimension and the number of beach maps.

        Args:
            file_name (str): The name of the file to be created.
            coeffs (List): The coefficients to be used in the LaTeX formula.
        """
        offset_a = 3 if int(2 * self.k) % 2 == 0 else 4
        offset_c = 5 if int(2 * self.k) % 2 == 0 else 7

        def format_binomial(coeff, j):
            if coeff == 1:
                return f'\\binom{{n - 1}}{{{offset_a + j}}}'
            return f'{coeff} \\binom{{n - 1}}{{{offset_a + j}}}'

        latex_formula = '= ' + ' \\\\\n & \\,\\,\\,\\, + '.join(
            (f'\\Bigg(' if sum(1 for x in row if x != 0) > 1 else '') +
            ' + '.join(format_binomial(coeff, j) for j, coeff in enumerate(row) if coeff != 0) +
            (f'\\Bigg)' if sum(1 for x in row if x != 0) > 1 else '') +
            f' (q - 1)^{{{offset_c + i}}}'
            for i, row in enumerate(coeffs)
        )

        k = self.k if int(2 * self.k) % 2 == 1 else int(self.k)

        output_file = Path(f"{file_name}/{file_name}.tex")
        output_file.parent.mkdir(exist_ok=True, parents=True)
        with open(output_file, "w") as file:
            file.write(self.latex_preamble())
            file.write(f"\\begin{{align*}}\n & \\dim(Z_{{{self.n, k}}}({self.q})) - |\\Pi_{{{self.n, k}}}({self.q})|"
                       f" \\\\ \n & {latex_formula} \n \\end{{align*}}\n")
            file.write("\\end{document}")

    def get_latex_table(self, file_name: str, coeffs: List[List[int]]):
        """
        Generates a LaTeX table displaying the coefficients for the equation of the difference of the dimension and
        the number of beach maps.

        Args:
            file_name (str): The name of the file to be created.
            coeffs: The coefficients to be displayed in the LaTeX table.
        """
        num_cols = len(coeffs[0])
        offset_a = 3 if int(2 * self.k) % 2 == 0 else 4
        offset_c = 5 if int(2 * self.k) % 2 == 0 else 7

        # Header
        latex_table = (
                "\\begin{table}[htb]\n"
                f"\t\\caption[Difference coefficients for k = {int(self.k) if int(2 * self.k) % 2 == 0 else self.k}]"
                f"{{\n"
                f"\t\tTable of values $A_{{i, j}}({int(self.k) if int(2 * self.k) % 2 == 0 else self.k})$. \n"
                "\t}\n"
                "\t\\begin{center}\n"
                "\t\t\\renewcommand{\\arraystretch}{1.5}\n"
                f"\t\t\\begin{{tabular}}{{||l{'|r' * num_cols}||}}\n\t\t\t\\hline \\hline\n\t\t\t "
                f"$k = {int(self.k) if int(2 * self.k) % 2 == 0 else self.k}$ & "
                f"" + " & ".join(f"$\\binom{{n - 1}}{{{offset_a + j}}}$" for j in range(num_cols)) + " \\\\ \\hline\n")

        # Rows
        for i, row in enumerate(coeffs):
            exponent = offset_c + i
            row_values = " & ".join(str(val) for val in row)
            latex_table += f"\t\t\t$(q - 1)^{{{exponent}}}$ & {row_values} \\\\ \\hline\n"

        # Footer
        latex_table += (
            "\t\t\t\\hline \\hline\n"
            "\t\t\\end{tabular}\n"
            f"\t\\label{{difftable{int(self.k) if int(2 * self.k) % 2 == 0 else self.k}}}\n"
            "\t\\end{center}\n"
            "\\end{table}"
        )

        output_file = Path(f"{file_name}/{file_name}.tex")
        output_file.parent.mkdir(exist_ok=True, parents=True)
        with open(output_file, "w") as file:
            file.write(self.latex_preamble())
            file.write(latex_table)
            file.write("\\end{document}")
