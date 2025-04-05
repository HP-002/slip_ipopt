def parse_file(filename):
    out_filename = "parsed_" + filename
    with open(filename, 'r') as input_file, open(out_filename, 'w') as output_file:
        for line in input_file:
            parsed_equation = modify_equation(line)
            output_file.write(parsed_equation)
    return True

def modify_equation(equation):
    # TODO: Modify the equations
    return equation


def main():
    input_files = ['input.txt', 'input2.txt']

    for file in input_files:
        parse_file(file)
    

if __name__ == '__main__':
    main()
