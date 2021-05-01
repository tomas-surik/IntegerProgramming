"""
OPTIMISATION FUNCTIONS

This code is based on my solution of an optimisation problem I had to solve.

It is an integer programming problem. The goal is to select options that 
results in the lowest costs. Number of options that must be selected is set. 
All available options have a set of three costs assigned to them. 

Costs A and B are independent across all options. However, C costs can be 
compensated across selected options. 

"""

import numpy as np
import cvxpy as cp
import pandas as pd

def optimise_total_costs(A_quantity, B_quantity, C_change, A_cost, B_cost, C_cost, C_base, num_options):
    """
    This function optimises total costs.
    """
    # Define the decision variable array
    x = cp.Variable(10, boolean = True)
    # Define constraints on decision variable
    constraints = [0 <= x, x <= 1, sum(x) == num_options]
    # Objective function 
    total_costs = cp.sum(A_cost * cp.multiply(A_quantity, x)) + cp.sum(B_cost * cp.multiply(B_quantity, x)) + (C_cost *  - cp.minimum(C_base + cp.sum(cp.multiply(C_change, x)), 0))
    # Goal of optimisation
    objective = cp.Minimize(total_costs)
    # Solve the optimisation problem 
    problem = cp.Problem(objective, constraints)
    result = problem.solve(solver = cp.ECOS_BB)
    return np.abs(np.rint(x.value))

def evaluate_total_costs(x, A_quantity, B_quantity, C_change, A_cost, B_cost, C_cost, C_base):
    """
    This function calculates total costs given the array of decision variables.
    """
    x = np.rint(x)
    A_total_costs = np.sum(A_cost * A_quantity * x)
    B_total_costs = np.sum(B_cost * B_quantity * x)
    C_total_costs = C_cost *  - min(C_base + np.sum(C_change * x), 0)
    return A_total_costs + B_total_costs + C_total_costs

def num_options_sensitivity(A_quantity, B_quantity, C_change, A_cost, B_cost, C_cost, C_base):
    """
    This code calculates optimal costs for very possible number of option that must be selected.
    """
    output = []
    for num_options in np.arange(1, 10 + 1):
        x_optimum = optimise_total_costs(A_quantity, B_quantity, C_change, A_cost, B_cost, C_cost, C_base, num_options)
        value_optimum = evaluate_total_costs(x_optimum, A_quantity, B_quantity, C_change, A_cost, B_cost, C_cost, C_base)
        output.append((num_options, x_optimum, value_optimum))
    return output

def AB_sensitivity(A_quantity, B_quantity, C_change, C_cost, C_base, num_options):
    """
    This code identifies optimal solutions for range of A and B unit cost values.
    """
    itin = 0
    A_cost_record = []
    B_cost_record = []
    x_optimum_record = []
    x_optimum_index_record = []
    value_optimum_record = []
    # Base values of costs
    A_cost_center = 300
    B_cost_center = 500
    # Largest increase or decrease that will be tested
    p1 = 0.7
    p2 = 0.5
    #  Step by which will by value of cost looped over
    s1 = 5
    s2 = 5
    A_cost_range = np.arange(np.rint(A_cost_center * (1 - p1)), A_cost_center * (1 + p1) + s1, s1)
    B_cost_range = np.arange(B_cost_center * (1 - p2), B_cost_center * (1 + p2) + s2, s2)
    print(A_cost_range, B_cost_range)
    for A_cost in A_cost_range:
        for B_cost in B_cost_range:
            print(A_cost, B_cost)
            itin += 1
            x_optimum = optimise_total_costs(A_quantity, B_quantity, C_change, A_cost, B_cost, C_cost, C_base, num_options)
            value_optimum = evaluate_total_costs(x_optimum, A_quantity, B_quantity, C_change, A_cost, B_cost, C_cost, C_base)
            A_cost_record.append(A_cost)
            B_cost_record.append(B_cost)
            x_optimum_record.append(x_optimum)
            x_optimum_index_record.append((np.nonzero(x_optimum) + np.array([1,1,1])).tolist())
            value_optimum_record.append(value_optimum) 
            print(itin)
    return pd.DataFrame.from_dict({"A_cost": A_cost_record, "B_cost": B_cost_record, "x_optimum": x_optimum_record, "flight_cancel": x_optimum_index_record, "costs": value_optimum_record})

def C_sensitivity(A_quantity, B_quantity, C_change, A_cost, B_cost, num_options):
    itin = 0
    C_cost_record = []
    C_base_record = []
    x_optimum_record = []
    x_optimum_index_record = []
    value_optimum_record = []
    # Base values of costs
    C_cost_center = 2500
    # Largest increase or decrease that will be tested
    p1 = 0.50
    #  Step by which will by value of cost looped over
    s1 = 50
    C_cost_range = np.arange(C_cost_center * (1 - p1), C_cost_center * (1 + p1) + s1, s1)
    # This value has floor equal to zero and only integer values are possible
    C_base_range = np.arange(0, 5 + 1, 1)
    print(C_cost_range, C_base_range)
    for C_cost in C_cost_range:
        for C_base in C_base_range:
            print(C_cost, C_base)
            itin += 1
            x_optimum = optimise_total_costs(A_quantity, B_quantity, C_change, A_cost, B_cost, C_cost, C_base, num_options)
            value_optimum = evaluate_total_costs(x_optimum, A_quantity, B_quantity, C_change, A_cost, B_cost, C_cost, C_base)
            C_cost_record.append(C_cost)
            C_base_record.append(C_base)
            x_optimum_record.append(x_optimum)
            x_optimum_index_record.append((np.nonzero(x_optimum) + np.array([1,1,1])).tolist())
            value_optimum_record.append(value_optimum) 
            print(itin)
    return pd.DataFrame.from_dict({"C_cost": C_cost_record, "C_base": C_base_record, "x_optimum": x_optimum_record, "flight_cancel": x_optimum_index_record, "costs": value_optimum_record})