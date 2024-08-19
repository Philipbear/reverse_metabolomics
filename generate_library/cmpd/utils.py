from molmass import Formula
import re


def neutralize_formula(formula):
    """
    deal with the charge in the formula
    such as C5H5N+ -> C5H4N
    """
    # Split the formula into the chemical part and the charge part
    match = re.match(r'([A-Za-z0-9]+)([-+]?\d*)', formula)
    if not match:
        return "Invalid formula"

    chemical, charge = match.groups()

    # If there's no charge, return the original formula
    if not charge:
        return chemical

    # Convert charge to integer
    if charge == '-':
        charge = -1
    elif charge == '+':
        charge = 1
    else:
        charge = int(charge) if charge else 0

    # Parse the chemical formula
    elements = re.findall(r'([A-Z][a-z]?)(\d*)', chemical)

    # Find H and its count
    h_index = next((i for i, (elem, _) in enumerate(elements) if elem == 'H'), None)
    if h_index is not None:
        h_count = int(elements[h_index][1]) if elements[h_index][1] else 1
    else:
        h_count = 0

    # Calculate new H count
    new_h_count = h_count - charge

    # Update or insert H in the elements list
    if new_h_count > 0:
        if h_index is not None:
            elements[h_index] = ('H', str(new_h_count) if new_h_count > 1 else '')
        else:
            elements.insert(1, ('H', str(new_h_count) if new_h_count > 1 else ''))
    elif h_index is not None:
        elements.pop(h_index)

    # Reconstruct the formula
    new_formula = ''.join(elem + count for elem, count in elements)

    return new_formula


def calc_exact_mass(formula):
    """
    Calculate the exact mass for a given formula string
    """
    try:
        f = Formula(formula)
        return f.monoisotopic_mass
    except:
        return None

