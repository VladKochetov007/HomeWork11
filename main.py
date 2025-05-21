import homework_module as hm
import math 

def run_demonstrations():
    output_lines = []

    # Task A
    output_lines.append("--- Task A: x_k = x^(2k) / (2k)! ---")
    test_cases_a = [(2.0, 3), (1.0, 5), (0.0, 4), (3.0, 0), (1.5, 4)]
    output_lines.append("Using generator (first 5 terms for x=2.0):")
    gen_a = hm.generate_x_k_terms(x=2.0)
    for i in range(5):
        try:
            output_lines.append(f"  x=2.0, k={i}: {next(gen_a):.8f}")
        except Exception as e:
            output_lines.append(f"  x=2.0, k={i}: Error - {e}")
            break


    output_lines.append("\nUsing loop function:")
    for x_val, k_val in test_cases_a:
        try:
            res = hm.calculate_x_k_loop(x_val, k_val)
            output_lines.append(f"  x={x_val}, k={k_val}: Result={res:.8f}")
        except Exception as e:
            output_lines.append(f"  x={x_val}, k={k_val}: Error - {e}")
    output_lines.append("-" * 30)

    # Task B
    output_lines.append("\n--- Task B: P_n = product_{i=1 to n} (1 + 1/i^2) ---")
    test_cases_b = [0, 1, 2, 3, 5, 10]
    output_lines.append("Using generator (P_0 to P_4):")
    gen_b = hm.generate_P_n_terms()
    for i in range(5): 
        try:
            output_lines.append(f"  n={i}: {next(gen_b):.8f}")
        except Exception as e:
            output_lines.append(f"  n={i}: Error - {e}")
            break

    output_lines.append("\nUsing loop function:")
    for n_val in test_cases_b:
        try:
            res = hm.calculate_P_n_loop(n_val)
            output_lines.append(f"  n={n_val}: Result={res:.8f}")
        except Exception as e:
            output_lines.append(f"  n={n_val}: Error - {e}")
    output_lines.append("-" * 30)

    # Task C
    output_lines.append("\n--- Task C: Determinant D_n ---")
    test_cases_c = [(1.0, 1.0, 3), (2.0, 3.0, 4), (1.0, 0.0, 5), (1.0,2.0,0), (1.0,2.0,1)]
    output_lines.append("Using generator (a=1, b=1, D_0 to D_4):")
    gen_c = hm.generate_D_n_terms(a=1.0, b=1.0)
    for i in range(5): 
        try:
            output_lines.append(f"  a=1,b=1, n={i}: {next(gen_c):.8f}")
        except Exception as e:
            output_lines.append(f"  a=1,b=1, n={i}: Error - {e}")
            break

    output_lines.append("\nUsing loop function:")
    for a_val, b_val, n_val in test_cases_c:
        try:
            res = hm.calculate_D_n_loop(a_val, b_val, n_val)
            output_lines.append(f"  a={a_val},b={b_val}, n={n_val}: Result={res:.8f}")
        except Exception as e:
            output_lines.append(f"  a={a_val},b={b_val}, n={n_val}: Error - {e}")
    output_lines.append("-" * 30)

    # Task D
    output_lines.append("\n--- Task D: S_n = sum_{k=1 to n} (a_k / 2^k) ---")
    test_cases_d = [0, 1, 2, 3, 4, 5, 10]
    output_lines.append("Using generator (partial sums S_1 to S_5):")
    gen_d = hm.generate_S_n_partial_sums()
    for i in range(1, 6): 
        try:
            output_lines.append(f"  n={i}: {next(gen_d):.8f}")
        except Exception as e:
            output_lines.append(f"  n={i}: Error - {e}")
            break
            
    output_lines.append("\nUsing loop function:")
    for n_val in test_cases_d:
        try:
            res = hm.calculate_S_n_loop(n_val)
            output_lines.append(f"  n={n_val}: Result={res:.8f}")
        except Exception as e:
            output_lines.append(f"  n={n_val}: Error - {e}")
    output_lines.append("-" * 30)
    
    # Task E
    output_lines.append("\n--- Task E: y = 2 * sum (x^(2m+1) / (2m+1)) ---")
    test_cases_e = [
        (0.5, 1e-6), (0.0, 1e-5), (-0.5, 1e-4), (0.9, 1e-7), (0.99, 1e-8)
    ]
    output_lines.append("Using generator (first 5 terms for x=0.5):")
    gen_e = hm.generate_taylor_sequence_terms_e(x=0.5)
    current_sum_e = 0.0
    for i in range(5):
        try:
            term = next(gen_e)
            current_sum_e += term
            output_lines.append(f"  x=0.5, term_{i}: {term:.8f}, current_sum: {current_sum_e:.8f}")
        except Exception as e:
            output_lines.append(f"  x=0.5, term_{i}: Error - {e}")
            break

    output_lines.append("\nUsing loop function:")
    for x_val, eps_val in test_cases_e:
        try:
            calc_y, math_y, terms = hm.calculate_y_e_loop(x_val, eps_val)
            output_lines.append(
                f"  x={x_val}, eps={eps_val}: Calculated y={calc_y:.8f} (in {terms} terms), Math lib y={math_y:.8f}"
            )
        except Exception as e:
            output_lines.append(f"  x={x_val}, eps={eps_val}: Error - {e}")
    output_lines.append("-" * 30)

    # Write to output.txt
    output_file_path = "output.txt"
    with open(output_file_path, "w", encoding="utf-8") as f:
        for line in output_lines:
            print(line) 
            f.write(line + "\n")
    
    print(f"\nDemonstration finished. Results are in {output_file_path} and printed above.")

if __name__ == "__main__":
    run_demonstrations()