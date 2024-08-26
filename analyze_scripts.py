import ast
import os
import argparse
import csv

def analyze_script(script_path: str) -> list[dict[str, str]]:
    """Analyze a Python script to extract information about its functions and dependencies."""
    with open(script_path, 'r') as file:
        tree = ast.parse(file.read(), script_path)
    
    functions = []
    imported_modules = set()

    for node in ast.walk(tree):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            for alias in node.names:
                imported_modules.add(alias.name.split('.')[0])
        elif isinstance(node, ast.FunctionDef):
            function_name = node.name
            args = [arg.arg for arg in node.args.args]
            returns = 'None' if node.returns is None else ast.dump(node.returns)
            dependencies = set()

            for subnode in ast.walk(node):
                if isinstance(subnode, ast.Call):
                    if isinstance(subnode.func, ast.Name) and subnode.func.id in imported_modules:
                        dependencies.add(subnode.func.id)
                    elif isinstance(subnode.func, ast.Attribute):
                        if isinstance(subnode.func.value, ast.Name) and subnode.func.value.id in imported_modules:
                            dependencies.add(subnode.func.value.id)

            functions.append({
                'Function Name': function_name,
                'Args': ', '.join(args),
                'Returns': returns,
                'Dependencies': ', '.join(dependencies)
            })

    return functions

def generate_report(directory: str) -> list[dict[str, str]]:
    """Generate a report of functions from all Python scripts in the specified directory."""
    report = []
    for filename in os.listdir(directory):
        if filename.endswith('.py'):
            script_path = os.path.join(directory, filename)
            functions = analyze_script(script_path)
            for func in functions:
                report.append({
                    'Script Name': filename,
                    'Function Name': func['Function Name'],
                    'Args': func['Args'],
                    'Returns': func['Returns'],
                    'Dependencies': func['Dependencies']
                })
    return report

def truncate_string(string: str, width: int) -> str:
    """Truncate a string if it exceeds the specified width, adding '...'."""
    return string if len(string) <= width else string[:width-3] + '...'

def calculate_column_widths(report: list[dict[str, str]]) -> dict[str, int]:
    """Calculate the maximum width needed for each column based on the report data."""
    headers = ['Script Name', 'Function Name', 'Args', 'Returns', 'Dependencies']
    col_widths = {header: len(header) for header in headers}
    for row in report:
        for key, value in row.items():
            col_widths[key] = max(col_widths[key], len(value))
    return col_widths

def print_report(report: list[dict[str, str]]):
    col_widths = {
        'Script Name': 20,
        'Function Name': 25,
        'Args': 40,
        'Returns': 100,
        'Dependencies': 15
    }

    # Header
    header = (
        f"| {truncate_string('Script Name', col_widths['Script Name'])} "
        f"| {truncate_string('Function Name', col_widths['Function Name'])} "
        f"| {truncate_string('Args', col_widths['Args'])} "
        f"| {truncate_string('Returns', col_widths['Returns'])} "
        f"| {truncate_string('Dependencies', col_widths['Dependencies'])} |"
    )
    print(header)
    print("|" + "|".join(['-' * col_widths[col] for col in col_widths]) + "|")  # Divider

    # Rows
    for row in report:
        print(
            f"| {truncate_string(row['Script Name'], col_widths['Script Name'])} "
            f"| {truncate_string(row['Function Name'], col_widths['Function Name'])} "
            f"| {truncate_string(row['Args'], col_widths['Args'])} "
            f"| {truncate_string(row['Returns'], col_widths['Returns'])} "
            f"| {truncate_string(row['Dependencies'], col_widths['Dependencies'])} |"
        )

def save_report_to_csv(report: list[dict[str, str]], csv_path: str):
    """Save the report to a CSV file."""
    with open(csv_path, 'w', newline='') as csvfile:
        fieldnames = ['Script Name', 'Function Name', 'Args', 'Returns', 'Dependencies']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for row in report:
            writer.writerow(row)

def save_report_to_markdown(report: list[dict[str, str]], md_path: str):
    """Save the report to a Markdown file."""
    if not report:
        print("No functions found in the given directory.")
        return
    
    col_widths = calculate_column_widths(report)
    
    with open(md_path, 'w') as mdfile:
        # Header
        header = (
            f"| {truncate_string('Script Name', col_widths['Script Name'])} "
            f"| {truncate_string('Function Name', col_widths['Function Name'])} "
            f"| {truncate_string('Args', col_widths['Args'])} "
            f"| {truncate_string('Returns', col_widths['Returns'])} "
            f"| {truncate_string('Dependencies', col_widths['Dependencies'])} |"
        )
        separator = f"|{'-' * col_widths['Script Name']}|{'-' * col_widths['Function Name']}|{'-' * col_widths['Args']}|{'-' * col_widths['Returns']}|{'-' * col_widths['Dependencies']}|"
        mdfile.write(header + "\n" + separator + "\n")

        # Rows
        for row in report:
            mdfile.write(
                f"| {truncate_string(row['Script Name'], col_widths['Script Name'])} "
                f"| {truncate_string(row['Function Name'], col_widths['Function Name'])} "
                f"| {truncate_string(row['Args'], col_widths['Args'])} "
                f"| {truncate_string(row['Returns'], col_widths['Returns'])} "
                f"| {truncate_string(row['Dependencies'], col_widths['Dependencies'])} |\n"
            )

def main():
    """Main function to parse arguments and generate reports."""
    parser = argparse.ArgumentParser(description="Generate a report of functions in Python scripts.")
    parser.add_argument('-d', '--directory', required=True, help='Path to the input files.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output CSV file.')
    parser.add_argument('-m', '--markdown', required=False, help='Path to the output Markdown file.')
    args = parser.parse_args()

    if not os.path.isdir(args.directory):
        print(f"Error: The directory {args.directory} does not exist.")
        return

    report = generate_report(args.directory)
    if not report:
        print("No Python scripts with functions found in the given directory.")
        return

    print_report(report)
    save_report_to_csv(report, args.output)

    if args.markdown:
        save_report_to_markdown(report, args.markdown)

if __name__ == '__main__':
    main()
