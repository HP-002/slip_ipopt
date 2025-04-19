import re

def modify_x_indices(line):
    # This regex matches x[n] where n is a number
    def update_index(match):
        index = match.group(1)
        if index == "0":
            return f"x[{index}]"
        else:
            return f"x[{index}+i]"

    return re.sub(r"x\[(\d+)\]", update_index, line)

def add_i(inFile, outFile):
    # Fix the indices for variable 'x'
    with open(inFile, 'r') as input_file, open(outFile, 'w') as output_file:
        for line in input_file:
            parsed_equation = modify_x_indices(line.rstrip('\n') + '\n')
            output_file.write(parsed_equation)
    return True

if __name__ == '__main__':
    input_file = 'Equations/equations_without_i.txt'
    output_file = 'Equations/equations_with_i.txt'
    add_i(input_file, output_file)
