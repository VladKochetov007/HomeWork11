import math
from typing import Generator

# Task A
# Recurrence: x_0 = 1, x_k = x_{k-1} * x^2 / ((2k)(2k-1)) for k >= 1

def generate_x_k_terms(x: float) -> Generator[float, None, None]:
    """
    Generates terms of the sequence x_k = x^(2k) / (2k)!.
    Recurrence: x_0 = 1
                x_k = x_{k-1} * x^2 / ((2k)(2k-1)) for k >= 1
    """
    if not isinstance(x, (int, float)):
        raise TypeError("Input x must be a number.")

    term = 1.0  # x_0
    yield term
    k_val = 1
    while True:
        try:
            denominator = (2 * k_val) * (2 * k_val - 1)
            if denominator == 0:
                term = 0.0 if x == 0 else float('inf') # Or handle as error
            else:
                term *= (x * x) / denominator
        except OverflowError:
            term = float('inf')
        except ZeroDivisionError: # Should be caught by denominator check, but as safeguard
            term = float('inf') if term > 0 else (float('-inf') if term < 0 else float('nan'))


        yield term
        k_val += 1


def calculate_x_k_loop(x: float, k_target: int) -> float:
    """
    Calculates x_k = x^(2k) / (2k)! for a specific k_target using a loop.
    Uses the recurrence: x_0 = 1, x_k = x_{k-1} * x^2 / ((2k)(2k-1)).
    """
    if not isinstance(x, (int, float)):
        raise TypeError("Input x must be a number.")
    if not isinstance(k_target, int) or k_target < 0:
        raise ValueError("k_target must be a non-negative integer.")

    if k_target == 0:
        return 1.0

    current_x_k = 1.0 # x_0
    for k_val in range(1, k_target + 1):
        denominator = (2 * k_val) * (2 * k_val - 1)
        if denominator == 0: 
             return float('nan') 
        if x == 0.0: 
            current_x_k = 0.0 
            break 
        current_x_k *= (x * x) / denominator
    return current_x_k

# Task B
# Recurrence: P_0 = 1, P_n = P_{n-1} * (1 + 1/n^2) for n >= 1

def generate_P_n_terms() -> Generator[float, None, None]:
    """
    Generates terms of the product P_n = product_{i=1 to n} (1 + 1/i^2).
    Yields P_0, P_1, P_2, ...
    Recurrence: P_0 = 1
                P_n = P_{n-1} * (1 + 1/n^2) for n >= 1
    """
    current_P_n = 1.0  # P_0
    yield current_P_n
    n_val = 1
    while True:
        current_P_n *= (1 + 1 / (n_val * n_val))
        yield current_P_n
        n_val += 1

def calculate_P_n_loop(n_target: int) -> float:
    """
    Calculates P_n = product_{i=1 to n} (1 + 1/i^2) for a specific n_target using a loop.
    Uses recurrence: P_0 = 1 (convention for empty product or starting point).
    The product in formula is from i=1. P_n = (1+1/1^2)...(1+1/n^2).
    So P_0 is 1 if we define it as the base for P_1 = P_0 * (1+1/1^2).
    If n_target is 0, result is 1 (empty product).
    """
    if not isinstance(n_target, int) or n_target < 0:
        raise ValueError("n_target must be a non-negative integer.")

    current_P_n = 1.0 
    for i_val in range(1, n_target + 1):
        current_P_n *= (1 + 1 / (i_val * i_val))
    return current_P_n


# Task C
# Recurrence: D_0 = 1, D_1 = a+b, D_n = (a+b)D_{n-1} - ab*D_{n-2} for n >= 2

def generate_D_n_terms(a: float, b: float) -> Generator[float, None, None]:
    """
    Generates terms of the determinant sequence D_n.
    Recurrence: D_0 = 1
                D_1 = a+b
                D_n = (a+b)D_{n-1} - ab*D_{n-2} for n >= 2
    """
    if not all(isinstance(val, (int, float)) for val in [a,b]):
        raise TypeError("Inputs a, b must be numbers.")

    D_prev2 = 1.0  # D_0
    yield D_prev2
    
    D_prev1 = a + b  # D_1
    yield D_prev1
    
    while True:
        D_current = (a + b) * D_prev1 - (a * b) * D_prev2
        yield D_current
        D_prev2, D_prev1 = D_prev1, D_current


def calculate_D_n_loop(a: float, b: float, n_target: int) -> float:
    """
    Calculates D_n for a specific n_target using a loop.
    Uses recurrence: D_0 = 1, D_1 = a+b, D_n = (a+b)D_{n-1} - ab*D_{n-2}.
    """
    if not all(isinstance(val, (int, float)) for val in [a,b]):
        raise TypeError("Inputs a, b must be numbers.")
    if not isinstance(n_target, int) or n_target < 0:
        raise ValueError("n_target must be a non-negative integer.")

    if n_target == 0:
        return 1.0
    if n_target == 1:
        return a + b

    D_prev2 = 1.0  # D_0
    D_prev1 = a + b  # D_1
    
    D_current = D_prev1 
    for _ in range(2, n_target + 1): 
        D_current = (a + b) * D_prev1 - (a * b) * D_prev2
        D_prev2, D_prev1 = D_prev1, D_current
    return D_current


# Task D
# S_n = sum_{k=1 to n} (a_k / 2^k)
# a_1=a_2=a_3=1, a_k = a_{k-1} + a_{k-3} for k >= 4

def _generate_a_k_task_d() -> Generator[int, None, None]:
    """
    Helper generator for sequence a_k for Task D.
    a_1=1, a_2=1, a_3=1
    a_k = a_{k-1} + a_{k-3} for k >= 4
    """
    a1, a2, a3 = 1, 1, 1
    yield a1
    yield a2
    yield a3
    
    q = [a1, a2, a3] # Stores [a_{k-3}, a_{k-2}, a_{k-1}]
    while True:
        a_k = q[2] + q[0] 
        yield a_k
        q.pop(0)
        q.append(a_k)

def generate_S_n_partial_sums() -> Generator[float, None, None]:
    """
    Generates partial sums S_n = sum_{k=1 to n} (a_k / 2^k).
    Yields S_1, S_2, S_3, ...
    Recurrence for S_n: S_0 = 0 (implicitly)
                        S_n = S_{n-1} + a_n / 2^n for n >= 1
    Uses helper _generate_a_k_task_d for a_k values.
    """
    a_k_generator = _generate_a_k_task_d()
    current_S_n = 0.0
    k_val = 1
    while True:
        try:
            a_k = next(a_k_generator)
        except StopIteration: 
            break
        term_k = a_k / (2**k_val)
        current_S_n += term_k
        yield current_S_n
        k_val += 1

def calculate_S_n_loop(n_target: int) -> float:
    """
    Calculates S_n = sum_{k=1 to n} (a_k / 2^k) for a specific n_target using a loop.
    a_k sequence: a_1=a_2=a_3=1, a_k = a_{k-1} + a_{k-3} for k >= 4.
    """
    if not isinstance(n_target, int) or n_target < 0:
        raise ValueError("n_target must be a non-negative integer.")

    if n_target == 0:
        return 0.0

    a_values = [0] * (n_target + 1) 

    current_S_n = 0.0
    for k in range(1, n_target + 1):
        if k == 1 or k == 2 or k == 3:
            a_values[k] = 1
        else: 
            a_values[k] = a_values[k-1] + a_values[k-3]
        
        current_S_n += a_values[k] / (2**k)
        
    return current_S_n


# Task E
# y = 2 * sum_{m=0 to inf} (x^(2m+1) / (2m+1)) with precision epsilon
# Recurrence for terms u_m = x^(2m+1)/(2m+1) in the sum (before multiplying by 2):
# u_0 = x
# u_m = u_{m-1} * x^2 * (2m-1)/(2m+1) for m >= 1

def generate_taylor_sequence_terms_e(x: float) -> Generator[float, None, None]:
    """
    Generates terms 2 * u_m = 2 * x^(2m+1) / (2m+1) for the Taylor series in Task E.
    y = sum of these terms.
    Recurrence for u_m: u_0 = x
                        u_m = u_{m-1} * x^2 * (2m-1)/(2m+1) for m >= 1
    """
    if not isinstance(x, (int, float)):
        raise TypeError("Input x must be a number.")
    if not (-1 < x < 1):
        raise ValueError("|x| must be strictly less than 1 for convergence.")

    u_m_component = x 
    yield 2 * u_m_component 

    m_val = 1 
    while True:
        factor = (x * x * (2 * m_val - 1)) / (2 * m_val + 1)
        if x == 0:
            u_m_component = 0.0
        elif abs(2*m_val + 1) < 1e-9: # Denominator is zero
             u_m_component = float('inf')
        else:
            u_m_component *= factor
        yield 2 * u_m_component
        m_val += 1


def calculate_y_e_loop(x: float, epsilon: float) -> tuple[float, float, int]:
    """
    Calculates y = 2 * sum_{m=0 to inf} (x^(2m+1) / (2m+1)) with precision epsilon.
    The sum stops when the absolute value of the current term to be added is less than epsilon.
    Returns (calculated_sum, math_library_value, number_of_terms_summed).
    """
    if not isinstance(x, (int, float)):
        raise TypeError("Input x must be a number.")
    if not (-1 < x < 1):
        raise ValueError("|x| must be strictly less than 1 for convergence.")
    if not isinstance(epsilon, float) or epsilon <= 0:
        raise ValueError("Epsilon must be a positive float.")

    current_y_sum = 0.0
    num_terms_summed = 0
    
    current_term_val = 2 * x # Term for m=0
    m = 0 
    max_iterations = 10000 

    while abs(current_term_val) >= epsilon and num_terms_summed < max_iterations:
        current_y_sum += current_term_val
        num_terms_summed += 1
        m += 1 

        if x == 0.0: 
            current_term_val = 0.0 
        else:
            numerator_factor = (x * x * (2 * m - 1))
            denominator_factor = (2 * m + 1)
            if abs(denominator_factor) < 1e-9: # Avoid division by zero
                current_term_val = float('inf') if current_term_val > 0 else float('-inf') # Or handle error
            else:
                current_term_val *= numerator_factor / denominator_factor
    
    if num_terms_summed == max_iterations and abs(current_term_val) >= epsilon:
        print(f"Warning (Task E): Max iterations ({max_iterations}) reached for x={x}, eps={epsilon}. Result may be inaccurate.")

    math_lib_val = float('nan')
    if abs(1.0 - x) < 1e-12: 
        math_lib_val = float('inf')
    elif abs(1.0 + x) < 1e-12:
        math_lib_val = float('-inf')
    else:
        try:
            math_lib_val = math.log((1 + x) / (1 - x))
        except (ValueError, ZeroDivisionError): 
             if x >= 1: math_lib_val = float('inf')
             elif x <= -1: math_lib_val = float('-inf')

    return current_y_sum, math_lib_val, num_terms_summed