import re

def parse_file(filename):
    out_filename = filename.rstrip('.txt') + '_parsed.txt'
    with open(filename, 'r') as input_file, open(out_filename, 'w') as output_file:
        for line in input_file:
            parsed_equation = parse_equation(line.rstrip('\n') + '\n')
            output_file.write(parsed_equation)
    return True

def parse_equation(equation):
    # Fix the indices for variable 'x'
    pattern = r'\bx(\d{1,2})\b'
    def replace(match):
        num = int(match.group(1)) - 1  # Subtract 1
        return f'x[{num}]'
    parsed_equation = re.sub(pattern, replace, equation)

    # TODO: Remove all instances of "Subscript()"
    parsed_equation = re.sub(r'Subscript\(\s*([A-Za-z0-9_]+)\s*,\s*([A-Za-z0-9_]+)\s*\)', r'\1\2', parsed_equation)

    # TODO: Remove all instances of "Power()"

    
    return parsed_equation


def main():
    input_files = ['Equations/Tests/test_equations_1.txt', 'Equations/Tests/test_equations_2.txt']

    for file in input_files:
        parse_file(file)
    

if __name__ == '__main__':
    main()
