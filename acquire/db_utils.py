from collections import Counter
import sqlite3
import random
from itertools import combinations, chain

DB = '../db/experiment_db.sqlite'

class Chemical:
    '''
    Object to process chemical component types; takes in chemical formulas, and stores dictionary, 'elements,'
    with counts of each element.
    '''
    ELEMENTS = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg'
    , 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
    'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr','Rb', 'Sr',
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb',
    'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
    'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir',
    'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm',
    'Bk','Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
     'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
    def __init__(self, formula=' ', notes='', molar_mass=0.0, price=0.0, is_salt = False):
        '''Create a new Chemical object. Each chemical object is uniquely identified by its evaluated formula.

        Keyword arguments:
        self --
        formula -- string of a given 'Chemical,' or electrolyte component. This can handle strings of the form Ca(ClO4)2, and will reorganize this to look like O8Cl2Ca, and store elements and their amounts in a dictionary format.
        notes -- some generic string that can be stored in the database about a given component. Typically this is used to mark if a given component is purchased in hydrate form, or to note the colloquial name.
        molar_mass -- molar mass is a float that is grams/mol
        price -- dollars/g if in grams; dollars/mL if in mL. This should be found assuming the maximum volume purchasable on sigma aldrich.
        is_salt -- boolean that denotes whether a component is a salt or a solvent.
        '''
        self.elements = self.parse_formula(formula)
        self.notes = notes
        self.molar_mass = molar_mass
        self.price = price # PRICE IS IN TERMS OF $/ML AND $/G
        self.is_salt = is_salt

    def parse_formula(self, formula): #recursive function to parse equivalent formulas
        '''Recursively parses through a formula string, accounting for parentheses and different orderings for elements.
        Returns a 'Counter' object, which is really just a dictionary with elements on the left, and amounts on the
        right.

        Keyword arguments:
        self -- 
        formula -- string of a given 'Chemical,' or electrolyte component.
        '''

        elements = Counter()
        i = 0
        while i < len(formula):
            if formula[i] == '(':
                count = 1
                for j in range(i + 1, len(formula)):
                    if formula[j] == '(':
                        count += 1
                    elif formula[j] == ')':
                        count -= 1
                        if count == 0:
                            break
                else:
                    raise TypeError(f'Unbalanced parentheses in formula {formula}')

                sub_elements = self.parse_formula(formula[i + 1:j])
                i = j + 1

                factor = ''
                while i < len(formula) and formula[i].isdigit():
                    factor += formula[i]
                    i += 1
                factor = int(factor) if factor else 1
                for element, quantity in sub_elements.items():
                    elements[element] += quantity * factor
            elif formula[i].isalpha():
                element = formula[i]
                i += 1
                while i < len(formula) and formula[i].islower():
                    element += formula[i]
                    i += 1
                if element not in self.ELEMENTS:
                    raise TypeError(f'Unknown element {element} in formula {formula}')
                quantity = ''
                while i < len(formula) and formula[i].isdigit():
                    quantity += formula[i]
                    i += 1
                quantity = int(quantity) if quantity else 1
                elements[element] += quantity
            else:
                raise TypeError(f'Invalid character {formula[i]} in formula {formula}')
        return elements

    def __eq__(self, other):
        '''Equivalence check which compares element dictionaries

        Keyword arguments:
        self -- 
        other -- the other Chemical object to be compared
        '''
        if isinstance(other, Chemical):
            return self.elements == other.elements
        return False

    def __str__(self):
        '''Outputs formula as a string with elements sorted by element number. Ex: Ca(ClO4)2 -> O8Cl2Ca

        Keyword arguments:
        self --
        '''
        sorted_elements = sorted(self.elements.items(), key=lambda x: self.ELEMENTS.index(x[0]))
        return ''.join(f'{element}{count}' if count > 1 else f'{element}' for element, count in sorted_elements)

def get_electrolyte_by_components(components: dict):
    '''Takes in dictionary of components and amounts, ex: {"formula1": 3.23, "formula2": .57} returns id of an electrolyte.

    Keyword arguments:
    components -- dictionary of components and amounts. Amounts should be a float, and components should be strings. Strings can be different to formulas listed in database, as a separate Chemical object is created for each formula to standardize formulas.
    '''
    conn = sqlite3.connect(DB, isolation_level='IMMEDIATE')
    try:
        c = conn.cursor()

        #creating long query string with all the requirements for amount and formula
        query_string = "SELECT e.id FROM electrolytes e WHERE e.id IN (SELECT ec.electrolyte_id FROM electrolyte_components ec WHERE "
        query_conditions = []
        query_values = []
        for formula, amount in components.items():
            chemical = Chemical(formula)
            query_conditions.append("(ec.component_id = (SELECT id FROM components WHERE formula = ?) AND ec.amount = ?)")
            query_values.extend([str(chemical), amount])
        query_string += " OR ".join(query_conditions) + " GROUP BY ec.electrolyte_id HAVING COUNT(ec.electrolyte_id) = ?)"
        query_values.append(len(components))

        c.execute(query_string, query_values)
        candidate_ids = [row[0] for row in c.fetchall()]

        #check that each candidate id doesn't have extra components
        
        for id in candidate_ids:
            c.execute("SELECT COUNT(*) FROM Electrolyte_Components WHERE Electrolyte_ID = ?", (id,))
            if c.fetchone()[0] != len(components):
                candidate_ids.remove(id)
        conn.commit()
    except sqlite3.Error as e:
        print(f"An error occurred: {e.args[0]}")
        conn.rollback()
    finally:
        conn.close()
    
    if(len(candidate_ids) > 1):
        raise ValueError(f"{len(candidate_ids)} total duplicate electrolytes found with components {str(components)}")
    return candidate_ids[0]

def add_electrolyte(components: dict,

                    conductivity: float,
                    conduct_uncert_bound: float,
                    concent_uncert_bound: float,

                    density: float = -1,
                    temperature: float = -1,
                    viscosity: float = -1,
                    v_window_low_bound: float = -1,
                    v_window_high_bound: float = -1,
                    surface_tension: float = -1
                    ):
    """adds a new electrolyte with certain properties--each electrolyte must be made up of componenets already stored in the database.

    Keyword arguments:
    components -- dict of chemical formula and amount; e.g. {str: float, ...}. Chemical objects will be created for each formula, so formulas do not have to conform exactly to those listed in component table.
    conductivity -- final conductivity value for given electrolyte--should NOT be resistance, conductance, or resistivity
    conduct_uncert_bound -- distance from uncertainty bounds for conductivity, i.e. \pm{conduct_uncert_bound}; -1 for no input.
    concent_uncert_bound -- distance from uncertainty bounds for concentration, i.e. \pm{concent_uncert_bound}; -1 for no input.
    density -- density, float, listed as -1 for no input.
    temperature -- temperature, float, listed as -1 for no input.
    viscosity -- viscocity, float, listed as -1 for no input.
    v_window_low_bound -- voltage window lower bound, float.
    v_window_high_bound -- voltage window upper bound, float.
    surface_tension -- surface tension, float, listed as -1 for no input.
    """
    if(check_electrolyte_exists(components)):
        raise ValueError(f'Electrolyte with formula {components} already exists')
    conn = sqlite3.connect(DB)
    try:
        c = conn.cursor()

        attr_dict = {
            'conductivity': conductivity,
            'conduct_uncert_bound': conduct_uncert_bound,
            'concent_uncert_bound': concent_uncert_bound,
            'density': density,
            'temperature': temperature,
            'viscosity': viscosity,
            'v_window_low_bound': v_window_low_bound,
            'v_window_high_bound': v_window_high_bound,
            'surface_tension': surface_tension,
        }

        c.execute(f'''INSERT INTO electrolytes {tuple(attr_dict.keys())} VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)'''
        , tuple(attr_dict.values()))

        conn.commit()

        #GET ID BACK FROM ELECTROLYTES TABLE BY CHECKING ALL ATTRIBUTES
        # Generate the parts of the WHERE clause
        clauses = [f"{attr} = ?" for attr in attr_dict.keys()]
        where_clause = " AND ".join(clauses)

        query = f"SELECT id FROM electrolytes WHERE {where_clause}"
        
        c.execute(query, tuple(attr_dict.values()))
        try:
            matched_ids = list(c.fetchall()[0])
        except:
            matched_ids = []

        c.execute('SELECT MAX(id) FROM electrolytes')
        electrolyte_id = c.fetchone()[0]

        if (electrolyte_id in matched_ids) and electrolyte_id != matched_ids:
            matched_ids.remove(electrolyte_id)
            print(f"Electrolytes with id(s): {matched_ids} have exactly identical attributes.")

        for formula, amount in components.items():
            chemical = Chemical(formula)
            c.execute("SELECT ID FROM components WHERE formula=?", (str(chemical),))

            component_id = c.fetchone()[0]
            print(electrolyte_id,component_id,amount)
            c.execute("INSERT INTO electrolyte_components (electrolyte_id, component_id, amount) VALUES (?, ?, ?)",
                    (electrolyte_id, component_id, amount))
            conn.commit()
    except sqlite3.Error as e:
        print(f"An error occurred: {e.args[0]}")
        conn.rollback()
    finally:
        conn.close()

def check_electrolyte_exists(components: dict):
    '''Checks if an electrolyte with matching components and amounts already exists in the database.

    Keyword arguments:
    components -- dict of chemical formula and amount; e.g. {str: float, ...}. Chemical objects will be created for each formula, so formulas do not have to conform exactly to those listed in component table.
    '''

    conn = sqlite3.connect(DB)
    try:
        c = conn.cursor()

        #long query string
        query_string = "SELECT e.id FROM electrolytes e WHERE e.id IN (SELECT ec.electrolyte_id FROM electrolyte_components ec WHERE "
        query_conditions = []
        query_values = []
        for formula, amount in components.items():
            chemical = Chemical(formula)
            query_conditions.append("(ec.component_id = (SELECT ID FROM components WHERE formula = ?) AND ec.amount = ?)")
            query_values.extend([str(chemical), amount])
        query_string += " OR ".join(query_conditions) + " GROUP BY ec.electrolyte_id HAVING COUNT(ec.electrolyte_id) = ?)"
        query_values.append(len(components))

        c.execute(query_string, query_values)
        candidate_ids = [row[0] for row in c.fetchall()]

        #check that each candidate id doesn't have extra components
        for id in candidate_ids:
            c.execute("SELECT COUNT(*) FROM Electrolyte_Components WHERE Electrolyte_ID = ?", (id,))
            if c.fetchone()[0] != len(components):
                candidate_ids.remove(id)
        conn.commit()
        # If we have at least one result, the electrolyte exists
        return len(candidate_ids) > 0
    except sqlite3.Error as e:
        print(f"An error occurred: {e.args[0]}")
        conn.rollback()
    finally:
        conn.close()

def get_components_by_id(electrolyte_id):
    '''takes in electrolyte id, and returns dictionary of components and amounts per component

    Keyword arguments:
    electrolyte_id -- int, electrolyte id
    '''
    conn = sqlite3.connect(DB)
    try:
        c = conn.cursor()

        query = """
        SELECT c.formula, ec.amount
        FROM electrolyte_components ec 
        JOIN components c ON ec.component_id = c.id
        WHERE ec.electrolyte_id = ?
        """

        c.execute(query, (electrolyte_id,))
        result = c.fetchall()
        # Return the components
        return {row[0]: row[1] for row in result}
    except sqlite3.Error as e:
        print(f"An error occurred: {e.args[0]}")
        conn.rollback()
    finally:
        conn.close()

def add_component_type(
    chemical: Chemical
):
    '''self-evident; only takes in Chemical class

    Keyword arguments:
    chemical -- the Chemical object in question.
    '''
    conn = sqlite3.connect(DB)
    try:
        c = conn.cursor()

        #check if chemical is already in database
        c.execute("SELECT * FROM components WHERE formula = ?", (chemical.__str__(),))
        conn.commit()
        rows = c.fetchall()

        if(len(rows) != 0):
            print(f'{str(len(rows))} entries with formula {chemical.__str__()} already in database')
        else:
            c.execute("INSERT INTO components (formula, notes, molar_mass, price, is_salt) VALUES (?,?,?,?,?)", (str(chemical),chemical.notes, chemical.molar_mass, chemical.price, chemical.is_salt))
            conn.commit()
    except sqlite3.Error as e:
        print(f"An error occurred: {e.args[0]}")
        conn.rollback()
    finally:
        conn.close()

def get_component_type(
    formula: str
):
    '''returns attributes about component from its formula.
    
    Keyword arguments:
    formula -- str, formula for component, does not have to be arranged; a Chemical object is created to reformat formula string when checking.
    '''
    formatted_formula = Chemical(formula).__str__()

    conn = sqlite3.connect(DB)
    try:
        c = conn.cursor()

        c.execute("SELECT * FROM components WHERE formula = ?", (Chemical(formula).__str__(),))
        conn.commit()
        rows = c.fetchall()

        if len(rows) > 1:
            print(f'Multiple components found for formula {formula}')
        elif len(rows) == 0:
            print(f'No components found for formula {formula}')
        else:
            #print(rows[0][1], rows[0][2], rows[0][3], rows[0][4])
            return Chemical(rows[0][1], rows[0][2], rows[0][3], rows[0][4])
    except sqlite3.Error as e:
        print(f"An error occurred: {e.args[0]}")
        conn.rollback()
    finally:
        conn.close()

def remove_component_type(
    formula: str
):
    '''removes component from table just from formula.

    Keyword arguments:
    formula -- str, formula for component, does not have to be arranged; a Chemical object is created to reformat formula string when checking.
    '''
    conn = sqlite3.connect(DB)
    try:
        c = conn.cursor()

        formatted = Chemical(formula).__str__()
        print(formatted)

        c.execute("SELECT * FROM components WHERE formula = ?", (formatted,))
        conn.commit()
        fetched = c.fetchall()
        if len(fetched) == 0:
            print(f'No components found for formula {formula}')
        else:
            c.execute("DELETE FROM components WHERE formula = ?", (formatted,))
            print("Deleted " + str(len(fetched)) + f" entries for formula {formula}.")
            conn.commit()
    except sqlite3.Error as e:
        print(f"An error occurred: {e.args[0]}")
        conn.rollback()
    finally:
        conn.close()

def remove_electrolyte_by_id(id:int):
    '''removes electrolyte entry, and corresponding electrolyte_components entries from electrolyte id

    Keyword arguments:
    id -- int id of electrolyte
    '''
    conn = sqlite3.connect(DB)
    try:
        c = conn.cursor()
        c.execute("DELETE FROM electrolytes WHERE id=?", (id,))
        c.execute("DELETE FROM electrolyte_components WHERE electrolyte_id=?", (id,))
        conn.commit()
    except sqlite3.Error as e:
        print(f"An error occurred: {e.args[0]}")
        conn.rollback()
    finally:
        conn.close()
def get_experiment_electrolytes(component_list):
    '''for a given experiment run, takes list of all considered components, and returns all electrolytes from database which are in the considered space (and subspaces). lists unincluded component amounts for all electrolytes in subspaces as 0.
    
    Keyword arguments:
    component_list -- list of formula strings.
    '''
    conn = sqlite3.connect(DB)
    try:
        c = conn.cursor()
        component_list = [str(Chemical(x)) for x in component_list] #convert to correct formatting
        #convert the list of components into a string format suitable for a SQL IN clause
        placeholders = ', '.join('?' for _ in component_list)
        c.execute(f'''SELECT formula FROM components WHERE formula IN ({placeholders})''', component_list)
        found_components = [row[0] for row in c.fetchall()]
        
        missing_components = set(component_list) - set(found_components)
        if missing_components:
            raise ValueError(f"The following components are not found in the database: {', '.join(missing_components)}")
        
        #get all electrolyte IDs that contain ANY of the provided components
        c.execute(f'''SELECT DISTINCT electrolyte_id
                    FROM electrolyte_components
                    WHERE component_id IN (SELECT id FROM components WHERE formula IN ({placeholders}))''', 
                component_list)
        possible_electrolyte_ids = [row[0] for row in c.fetchall()]
        
        #filter the above results to get only those electrolyte IDs that contain ONLY the provided components or a subset of them
        valid_electrolyte_ids = []

        c.execute(f'''SELECT id FROM components WHERE formula IN ({placeholders})''', component_list)
        component_ids = [row[0] for row in c.fetchall()]

        for eid in possible_electrolyte_ids:
            c.execute('''SELECT DISTINCT component_id
                        FROM electrolyte_components
                        WHERE electrolyte_id = ?''', (eid,))
            components_of_eid = [row[0] for row in c.fetchall()]

            #print(f"For Electrolyte ID {eid}, Components are: {components_of_eid}")
            #check if all the components of the electrolyte are within the provided list
            if all(comp_id in component_ids for comp_id in components_of_eid):
                #print(f"Electrolyte ID {eid} is valid")
                valid_electrolyte_ids.append(eid)
            #else:
                #print(f"Electrolyte ID {eid} is NOT valid")

        #fetch the details of the valid electrolytes and their components
        results = []
        for eid in valid_electrolyte_ids:
            c.execute('SELECT * FROM electrolytes WHERE id = ?', (eid,))
            electrolyte_data = c.fetchone()

            #for each component in the provided list, get its amount or default to 0
            amounts = []
            for comp in component_list:
                c.execute('''
                    SELECT amount 
                    FROM electrolyte_components 
                    JOIN components ON components.id = electrolyte_components.component_id
                    WHERE electrolyte_id = ? AND formula = ?''', (eid, comp))
                row = c.fetchone()
                if row:
                    amounts.append(row[0])
                else:
                    amounts.append(0)
            
            results.append({
                'electrolyte': electrolyte_data,
                'amounts': amounts
            })
        return results

    except sqlite3.Error as e:
        print(f"An error occurred: {e.args[0]}")
        conn.rollback()
    finally:
        conn.close()


def generate_test_electrolytes(num_electrolytes: int):
    '''generates dummy electrolytes with randomized conductivities; inputs these into database.
    
    Keyword arguments:
    num_electrolytes -- number of dummies to generate.
    '''
    def generate_random_amount(upper_limit=31.9):
        return random.uniform(0.1, upper_limit)

    def calculate_volume(components_dict, components_info_dict):
        return sum(amount for formula, amount in components_dict.items() if not components_info_dict[formula])
    conn = sqlite3.connect(DB)
    c = conn.cursor()

    # Get all component formulas and their salt status
    c.execute("SELECT formula, is_salt FROM components")
    components_info = {formula: is_salt for formula, is_salt in c.fetchall()}

    # Generate all possible combinations of at least two component formulas
    all_combinations = list(chain.from_iterable(combinations(components_info.keys(), r) for r in range(2, len(components_info)+1)))
    
    created_electrolytes = 0

    # Step 1 & 2: Use unique combinations
    for combo in all_combinations:
        if created_electrolytes >= num_electrolytes:
            break
        
        components_dict = {formula: generate_random_amount() for formula in combo}
        while calculate_volume(components_dict, components_info) > 32:
            components_dict = {formula: generate_random_amount() for formula in combo}

        add_electrolyte(components_dict, random.uniform(0,10), 0.1, 0.1)  # Using dummy values for conductivity, etc.
        created_electrolytes += 1

    # Step 3: Generate additional random electrolytes
    while created_electrolytes < num_electrolytes:
        num_components = random.randint(2, len(components_info))
        combo = random.sample(components_info.keys(), num_components)
        
        components_dict = {formula: generate_random_amount() for formula in combo}
        while calculate_volume(components_dict, components_info) > 32:
            components_dict = {formula: generate_random_amount() for formula in combo}

        add_electrolyte(components_dict, 1.0, 0.1, 0.1)  # Using dummy values for conductivity, etc.
        created_electrolytes += 1

    conn.close()