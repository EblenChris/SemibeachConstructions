from CentralizerAlgebra import *
from Semibeaches import *
from sympy import pprint
import ast
import sys


def safe_input(prompt):
    """Handles user input with an option to exit."""
    try:
        user_input = input(prompt).strip().lower()
        if user_input == "exit":
            print("Exiting...")
            sys.exit(0)
        return user_input
    except KeyboardInterrupt:
        print("\nExit requested. Exiting...")
        sys.exit(0)


def get_user_input(choice: str):
    """Prompts the user to input values for n, k, and q, with exit options."""
    if choice == "random_semibeach" or choice == "two_beaches":
        while True:
            try:
                n = int(safe_input("Enter a positive integer n > 1: "))
                if n > 1:
                    break
                print("Invalid input.")
            except ValueError:
                print("Invalid input.")

        while True:
            try:
                k = float(safe_input("Enter a positive integer or half-integer k: "))
                if k > 0 and 2 * k == int(2 * k):
                    break
                print("Invalid input.")
            except ValueError:
                print("Invalid input.")

        while True:
            try:
                q = int(safe_input("Enter a positive prime power q: "))
                if q > 0:
                    break
                print("Invalid input.")
            except ValueError:
                print("Invalid input.")
        return n, k, q
    elif choice == "difference_coefficients":
        while True:
            try:
                k = float(safe_input("Enter a positive integer or half-integer k: "))
                if k > 0 and 2 * k == int(2 * k):
                    break
                print("Invalid input.")
            except ValueError:
                print("Invalid input.")
        return int(2 * k), k, 2
    elif choice == "spreadsheet":
        while True:
            try:
                n_values = convert_input_to_list(safe_input("Enter a list, e.g. [1, 2, 3], of positive "
                                                            "integers n > 1: "), "spreadsheet")
                if len(n_values) > 0 and all(n > 1 for n in n_values):
                    break
                print("Invalid input.")
            except ValueError:
                print("Invalid input.")

        while True:
            try:
                k_values = convert_input_to_list(safe_input("Enter a list, e.g. "
                                                            "[1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6], of positive"
                                                            " integers or half-integers k: "), "spreadsheet")
                if len(k_values) > 0 and all(k > 0 and 2 * k == int(2 * k) for k in k_values):
                    break
                print("Invalid input.")
            except ValueError:
                print("Invalid input.")

        while True:
            try:
                q_values = convert_input_to_list(safe_input("Enter a list, e.g. [2], of positive prime "
                                                            "powers q: "), "spreadsheet")
                if len(q_values) > 0 and all(q > 0 for q in q_values):
                    break
                print("Invalid input.")
            except ValueError:
                print("Invalid input.")
        return n_values, k_values, q_values


def generate_two_beaches():
    n, k, q = get_user_input("two_beaches")
    while True:
        try:
            choice = safe_input("\nChoose an option: "
                                "\n1. Enter a labeled set partition"
                                "\n2. Generate a random labeled set partition"
                                "\n")
            if choice == '1':
                try:
                    set_partition_str = safe_input(f"Paste a valid set partition of {n} here: ")
                    set_partition = convert_input_to_list(set_partition_str, "specific_set_partition")
                    if (ArcDiagram(set_partition, n, q).is_set_partition() and
                            (n not in ArcDiagram(set_partition, n, q).rights or int(2 * k) % 2 == 0)):
                        break
                    print("Invalid input.")
                except ValueError:
                    print("Invalid input.")
            elif choice == "2":
                while choice == "2":
                    if int(2 * k) % 2 == 0:
                        set_partitions = SetPartitionGenerator(n, q).labeled_set_partitions()
                    else:
                        set_partitions = SetPartitionGenerator(n - 1, q).labeled_set_partitions()
                    super_chars = [sp for sp in set_partitions if len(sp) <= floor(k)]
                    set_partition = random.choice(super_chars)
                    while True:
                        try:
                            choice = safe_input(f"\nYour labeled set partition is {set_partition}."
                                                f"\n\nChoose an option."
                                                f"\n1. Generate two beaches with this final set partition."
                                                f"\n2. Generate a different labeled set partition."
                                                f"\n")
                            if choice in ['1', '2']:
                                break
                            else:
                                print("Invalid input.")
                        except ValueError:
                            print("Invalid input.")
                break
        except ValueError:
            print("Invalid input.")
    print("\nGenerating two beaches with given final set partition...")
    beach_1, beach_2 = (BeachGenerator(n, k, q).random_beach(set_partition),
                        BeachGenerator(n, k, q).random_beach(set_partition))
    file_name = "two_beaches"
    Semibeach(beach_1, q).tex_draw(file_name, beach_2)
    print(f"\nCheck the folder to view the {file_name}.tex file.")


def convert_input_to_list(input_str: str, choice: str):
    """Converts a string representation of a list of into a list object."""
    if choice == "specific_set_partition":
        try:
            result = ast.literal_eval(input_str)
            if isinstance(result, list) and all(
                isinstance(item, tuple) and len(item) == 2 and
                isinstance(item[0], tuple) and len(item[0]) == 2 and
                all(isinstance(x, int) for x in item[0])
                for item in result
            ):
                return result
            else:
                raise ValueError("Input does not match the required format.")
        except (ValueError, SyntaxError):
            raise ValueError("Invalid input format.")
    elif choice == "spreadsheet":
        try:
            result = ast.literal_eval(input_str)  # Safely evaluate input
            if isinstance(result, list) and all(isinstance(x, (int, float)) for x in result):
                return result
            else:
                raise ValueError("Input must be a list of numbers.")
        except (ValueError, SyntaxError):
            raise ValueError("Invalid input format.")


def generate_random_semibeach():
    """Generates a random semibeach along with its associated ."""
    n, k, q = get_user_input("random_semibeach")

    while True:
        want_beach = safe_input("Do you want a beach or a semibeach? Type 'beach' or 'semibeach': ")
        if want_beach in {'beach', 'semibeach'}:
            break
        print("Invalid input.")

    if want_beach == 'beach':
        while True:
            try:
                choice = safe_input("\nChoose an option: "
                                    "\n1. Enter a labeled set partition"
                                    "\n2. Generate a random labeled set partition"
                                    "\n")
                if choice == '1':
                    try:
                        set_partition_str = safe_input(f"Paste a valid set partition of {n} here: ")
                        set_partition = convert_input_to_list(set_partition_str, "specific_set_partition")
                        if (ArcDiagram(set_partition, n, q).is_set_partition() and
                                (n not in ArcDiagram(set_partition, n, q).rights or int(2 * k) % 2 == 0)):
                            break
                        print("Invalid input.")
                    except ValueError:
                        print("Invalid input.")
                elif choice == "2":
                    while choice == "2":
                        if int(2 * k) % 2 == 0:
                            set_partitions = SetPartitionGenerator(n, q).labeled_set_partitions()
                        else:
                            set_partitions = SetPartitionGenerator(n - 1, q).labeled_set_partitions()
                        super_chars = [sp for sp in set_partitions if len(sp) <= floor(k)]
                        set_partition = random.choice(super_chars)
                        while True:
                            try:
                                choice = safe_input(f"\nYour labeled set partition is {set_partition}."
                                                    f"\n\nChoose an option."
                                                    f"\n1. Generate a random beach with this final set partition."
                                                    f"\n2. Generate a different labeled set partition."
                                                    f"\n")
                                if choice in ['1', '2']:
                                    break
                                else:
                                    print("Invalid input.")
                            except ValueError:
                                print("Invalid input.")
                    break
            except ValueError:
                print("Invalid input.")
        print("\nGenerating random beach with given final set partition...")
        b, shifts = BeachGenerator(n, k, q).random_beach(set_partition)
    else:
        print("\nGenerating random semibeach...")
        b, shifts = SemibeachGenerator(n, k, q).random_semibeach(False)

    file_name = 'beach' if want_beach == 'beach' else 'semibeach'
    Semibeach((b, shifts), q).tex_draw(file_name, None)
    print(f"\nCheck the folder to view the {file_name}.tex file.")


def compute_difference_coefficients():
    """Computes and displays the coefficients for the difference formula."""
    n, k, q = get_user_input("difference_coefficients")
    print("\nComputing difference coefficients...")
    coeff_matrix, sheets_syntax, latex_syntax = CentralizerAlgebra(n, k, q).get_difference_coeffs()
    print("\nCoefficient Matrix:")
    pprint(coeff_matrix)
    print("\nGoogle Sheets Syntax:")
    print(sheets_syntax)
    print("\nLaTeX Syntax:")
    print(latex_syntax)
    print(f"\nCheck the folder for the centralizer_algebra_results.csv file.")


def get_spreadsheet_of_values():
    n_values, k_values, q_values = get_user_input("spreadsheet")
    n, k, q = n_values[0], k_values[0], q_values[0]
    extra_n, extra_k, extra_q = n_values[1:], k_values[1:], q_values[1:]
    CentralizerAlgebra(n, k, q).get_spreadsheet(extra_n, extra_k, extra_q)
    print(f"\nCheck the folder for the centralizer_algebra_results.csv file.")


def main():
    """Main function to run the desired computations interactively, with exit options."""
    while True:
        print("\nChoose an option. You may exit at any time by typing 'exit': ")
        print("1. Generate random semibeach.")
        print("2. Generate two beaches with the same final set partition.")
        print("3. Compute difference coefficients.")
        print("4. Get spreadsheet of values of dimension, number of beach maps, and their ratio for a range of Args.")

        choice = safe_input("Enter the number of your choice: ")
        if choice == '1':
            generate_random_semibeach()
        elif choice == '2':
            generate_two_beaches()
        elif choice == '3':
            compute_difference_coefficients()
        elif choice == '4':
            get_spreadsheet_of_values()
        elif choice == 'exit':
            print("Exiting.")
            break
        else:
            print("Invalid choice. Please select a valid option.")


if __name__ == "__main__":
    main()
