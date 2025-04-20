import re

def parse_file(inFile, outFile):
    with open(inFile, 'r') as input_file, open(outFile, 'w') as output_file:
        for line in input_file:
            parsed_equation = parse_equation(line.rstrip('\n') + '\n')
            output_file.write(parsed_equation)
    return True

def parse_equation(equation):
    # Fix the indices for variable 'x'
    def replace(match):
        num = int(match.group(1)) - 1  # Subtract 1
        return f'x[{num}]'
    
    equation = re.sub(r'\bx(\d{1,2})\b', replace, equation)

    # Remove all instances of "Subscript()"
    equation = re.sub(r'Subscript\(\s*([A-Za-z0-9_]+)\s*,\s*([A-Za-z0-9_]+)\s*\)', r'\1\2', equation)

    # Replace all instances of "Power(x,y)" with "pow(x,y)"
    equation = re.sub(r'\bPower\b', 'pow', equation)

    # Replace all instances of "Sin(x)" with "sin(x)"
    equation = re.sub(r'\bSin\b', 'sin', equation)

    # Replace all instances of "Cos(x)" with "cos(x)"
    equation = re.sub(r'\bCos\b', 'cos', equation)

    # Replace all constants to avoid variable name conflicts between IPOpt and equations
    equation = re.sub(r'\bg\b', 'G_', equation)     # g -> G_
    equation = re.sub(r'\bRa\b', 'Ra_', equation)   # Ra -> Ra_
    equation = re.sub(r'\bR\b', 'R_', equation)     # R -> R_
    equation = re.sub(r'\bJ\b', 'J_', equation)     # J -> J_
    equation = re.sub(r'\bb\b', 'B_', equation)     # b -> B_
    equation = re.sub(r'\bc\b', 'C_', equation)     # c -> C_
    equation = re.sub(r'\bm\b', 'M_', equation)     # m -> M_
    equation = re.sub(r'\bl0\b', 'L0_', equation)   # l0 -> L0_
    equation = re.sub(r'\bLa\b', 'La_', equation)   # La -> La_
    equation = re.sub(r'\bk0\b', 'K0_', equation)   # k0 -> K0_
    equation = re.sub(r'\bkb\b', 'Kb_', equation)   # kb -> Kb_
    equation = re.sub(r'\bkt\b', 'Kt_', equation)   # kt -> Kt_

    return equation


def main():
    # Test files
    # input_files = ['Equations/Tests/test_equations_1.txt', 'Equations/Tests/test_equations_2.txt', 'Equations/Tests/test_equations_3.txt', 'Equations/Tests/test_equations_4.txt']

    # Equations files for Col_var_v
    # input_files = ['Equations/Col_var_v/cost_var_v_text.txt', 'Equations/Col_var_v/costgrad_var_v_text.txt', 'Equations/Col_var_v/equality_var_v_text.txt', 'Equations/Col_var_v/inequality_var_v_text.txt', 'Equations/Col_var_v/eq_jacob_var_v_text.txt'];
    # output_files = ['Equations/Col_var_v/parsed/cost_var_v_parsed.txt', 'Equations/Col_var_v/parsed/costgrad_var_v_parsed.txt', 'Equations/Col_var_v/parsed/equality_var_v_parsed.txt', 'Equations/Col_var_v/parsed/inequality_var_v_parsed.txt', 'Equations/Col_var_v/parsed/eq_jacob_var_v_parsed.txt'];
    input_files = ['Equations/Col_var_v/obj_hessian_var_v_text.txt'];
    output_files = ['Equations/Col_var_v/parsed/obj_hessian_var_v_parsed.txt'];

    for inFile, outFile in zip(input_files, output_files):
        parse_file(inFile, outFile)
    

if __name__ == '__main__':
    main()
