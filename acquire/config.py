import json
import os
from typing import Tuple, List
from datetime import datetime as dt

from db_utils import Chemical

def generate_config_json(all_chemicals: Tuple[Chemical, Tuple[float,float]],
                        used_chemicals: List[Chemical],
                        max_vol: float,
                        min_vol: float,
                        initials = '',
                        output_filepath = '.'
                        ):
    '''
    bounds list is of form:
    {
    Chemical: tuple[lower_bound_float, upper_bound_float]
    }
    '''

    DISCRETE_RESOLUTION = .001

    domain = {}

    for chemical, bounds in all_chemicals:
        key = f"{str(chemical)}_grams" if chemical.is_salt else f"{str(chemical)}_mLs"
        domain[key] = {
            "name": key,
            "type": "discrete_numeric",
            "items": f"{bounds[0]}:{str(DISCRETE_RESOLUTION)}:{bounds[1]}"
        }

    now = dt.isoformat(dt.now()).replace(':','-')
    now = now[:now.index('.')]

    CONSTRAINT_FILENAME = 'constr_' + now + initials +'.py'
    EXP_NAME = 'exp_' + now + initials

    def generate_bool_str_unused_chemicals() -> str:
        unused_chemicals = [chemical for chemical, bounds in all_chemicals if chemical not in used_chemicals]
        return_str = ''

        constraint_strings = [f"{str(chemical)}_grams" if chemical.is_salt else f"{str(chemical)}_mLs" + '==0'
                              for chemical in unused_chemicals
                              ]
        return_str += constraint_strings.pop()
        for string in constraint_strings:
            return_str += ' and ' + string
        
        return '(' + return_str + ')'
    used_chemical_strings = [str(used) for used in used_chemicals]
    vol_args = [f"{str(chemical)}_mLs"
                    if
                        (not chemical.is_salt)
                        and str(chemical) in used_chemical_strings
                    else ''
                    for chemical in used_chemicals
                    ]
    min_vol_str = '(' + '+'.join(arg for arg in vol_args) + '>=' + str(min_vol) + ')'
    max_vol_str = '(' + '+'.join(arg for arg in vol_args) + '<=' + str(max_vol) + ')'
    unused_bool_str = generate_bool_str_unused_chemicals()

    full_constraints_str = ' and '.join([unused_bool_str, min_vol_str, max_vol_str])
    
    constraint_function_string = f'''
#Constraint file for experiment: {EXP_NAME}

def constraint(raw_point):
    {', '.join([f"{str(chemical)}_grams" if chemical.is_salt else f"{str(chemical)}_mLs" for chemical,bounds in all_chemicals])} = raw_point
    return {full_constraints_str}
'''
    #print(constraint_function_string)
    result = {
        "name": EXP_NAME,
        "domain": domain,
        "domain_constraints":{
            "constraint_1":{
                "name":"combined_constraints",
                "constraint":CONSTRAINT_FILENAME
            }
        }
    }
    with open(os.path.join(output_filepath, EXP_NAME + '.json'), 'w') as file:
        json.dump(result, file, indent=2)
    with open(os.path.join(output_filepath, CONSTRAINT_FILENAME), 'w') as file:
        file.write(constraint_function_string)
    
    return output_filepath +'/'+ EXP_NAME + '.json', output_filepath +'/'+ CONSTRAINT_FILENAME

sample_chemicals = [
    (Chemical("Na2SO4", is_salt=True), (0, 7.25)),
    (Chemical("NaNO3"), (0, 7.25)),
    (Chemical("LiNO3"), (0, 7.25)),
    (Chemical("Li2SO4"), (0, 7.25)),
    (Chemical("NaClO4"), (0, 7.25))
]
"""
used_chemicals = [Chemical("Na2SO4"),Chemical("NaNO3"),Chemical("LiNO3")]
generate_config_json(sample_chemicals, used_chemicals,32.0,0.0)
"""