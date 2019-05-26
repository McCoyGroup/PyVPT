"""Simple utilities that support the parsing part of pulling stuff out of Gaussian .log files
"""

import numpy as np, re

# we'll define a pile of Regex strings to use as components when designing matches
# for pulling stuff out of the blocks that we get from the GaussianReader

# wrapper patterns
grp_p = lambda p: r"("+p+r")" # capturing group
non_cap_p = lambda p: r"(?:"+p+r")" # non-capturing group
op_p = lambda p: r"("+p+r")?" # optional group
opnb_p = lambda p: r"(?:"+p+r")?" # optional non-binding group
rep_p = lambda p, n, m: r"("+p+"){"+str(n)+","+str(m)+"}"
repnb_p = lambda p, n, m: r"(?:"+p+"){"+str(n)+","+str(m)+"}"

sign_p = r"[\+\-]"
paren_p = r"\("+".*?"+"\)"
num_p = opnb_p(sign_p)+r"\d*\.\d+" # real number
int_p = opnb_p(sign_p)+r"\d+" # integer
posint_p = r"\d+" # only positive integer
ascii_p = "[a-zA-Z]"
name_p = ascii_p+"{1,2}" # atom name
ws_char_class = r"(?!\n)\s" # probably excessive... but w/e I'm not winning awards for speed here
ws_p = ws_char_class+"*" # whitespace
wsr_p = ws_char_class+"+" # real whitespace
cart_p = ws_p.join([ grp_p(num_p) ]*3) # cartesian coordinate
acart_p = "("+int_p+")"+ws_p+cart_p # atom coordinate as comes from a XYZ table

acart_p_c = re.compile(acart_p) # adds a little bit of a performance boost
cart_p_c = re.compile(cart_p)
cart_str_c = re.compile(ws_p.join([ num_p ]*3))
def pull_coords(txt, regex = cart_str_c, coord_dim = 3):
    """Pulls just the Cartesian coordinates out of the string

    :param txt: some string that has Cartesians in it somewhere
    :type txt:
    :param regex: the Cartesian matching regex; can be swapped out for a different matcher
    :type regex:
    :param coord_dim: the dimension of the pulled coordinates (used to reshape)
    :type coord_dim:
    :return:
    :rtype:
    """
    coords = re.findall(regex, txt)
    base_arr = np.array(coords, dtype=np.str).astype(dtype=np.float64)
    num_els = base_arr.shape[0]
    new_shape = (int(num_els/coord_dim), coord_dim) # for whatever weird reason this wasn't int...?
    new_arr = np.reshape(base_arr, new_shape)
    return new_arr

def pull_xyz(txt, num_atoms = None, regex = acart_p_c):
    """Pulls XYX-type coordinates out of a string

    :param txt: XYZ string
    :type txt: str
    :param num_atoms: number of atoms if known
    :type num_atoms: int
    :return: atom types and coords
    :rtype:
    """
    if num_atoms is None:
        num_cur = 15 # probably more than we'll actually need...
        atom_types = [None]*num_cur
        coord_array = np.zeros((num_cur, 3))
        for i, match in enumerate(re.finditer(regex, txt)):
            if i == num_cur:
                atom_types.extend([None]*(2*num_cur))
                coord_array = np.concatenate((coord_array, np.zeros((2*num_cur, 3))), axis=1)
            g = match.groups()
            atom_types[i] = g[:-3]
            coord_array[i] = np.array(g[-3:], dtype=np.str).astype(np.float64)

    else:
        atom_types = [None]*num_atoms
        coord_array =  np.zeros((num_atoms, 3), dtype=np.float64)
        parse_iter = re.finditer(regex, txt)
        for i in range(num_atoms):
            match = next(parse_iter)
            g = match.groups()
            atom_types[i] = g[:-3]
            coord_array[i] = np.array(g[-3:], dtype=np.str).astype(np.float64)
    atom_types = atom_types[:i+1]
    coord_array = coord_array[:i+1]

    return (atom_types, coord_array)

num_p_c = re.compile(num_p)
def pull_zmat_coords(txt, regex = num_p_c):
    '''Pulls only the numeric coordinate parts out of a Z-matrix (drops ordering or types)

    :param txt:
    :type txt:
    :param regex:
    :type regex:
    :return:
    :rtype:
    '''
    coords = re.findall(regex, txt)
    # print(txt)
    # print(coords)
    base_arr = np.array(coords, dtype=np.str).astype(dtype=np.float64)
    num_els = base_arr.shape[0]
    if num_els == 1:
        base_arr = np.concatenate((base_arr, np.zeros((2,), dtype=np.float64)))
        num_els = 3
    elif num_els == 3:
        base_arr = np.concatenate((base_arr, np.zeros((1,), dtype=np.float64)))
        base_arr = np.insert(base_arr, 1, np.zeros((2,), dtype=np.float64))
        num_els = 6
    else:
        base_arr = np.insert(base_arr, 3, np.zeros((1,), dtype=np.float64))
        base_arr = np.insert(base_arr, 1, np.zeros((2,), dtype=np.float64))
        num_els = num_els + 3

    # print(base_arr)
    coord_dim = 3
    new_shape = (int(num_els/coord_dim), coord_dim)
    new_arr = np.reshape(base_arr, new_shape)
    return new_arr

def process_zzzz(i, g, atom_types, index_array, coord_array, num_header=1):
    g_num = len(g)
    # print(g, g_num, num_header)
    if g_num == num_header+0:
        atom_types[i] = g
    elif g_num == num_header+2:
        atom_types[i] = g[:-2]
        # make atom refs to insert into array
        ref = np.array(g[-2:-1], dtype=np.str).astype(np.int8)
        ref = np.concatenate((ref, np.zeros((2,), dtype=np.int8)))
        coord = np.array(g[-1:], dtype=np.str).astype(np.float64)
        coord = np.concatenate((coord, np.zeros((2,), dtype=np.float64)))
        index_array[i-1] = ref
        coord_array[i-1] = coord
    elif g_num == num_header+4:
        atom_types[i] = g[:-4]
        # make atom refs to insert into array
        ref = np.array(g[-4::2], dtype=np.str).astype(np.int8)
        ref = np.concatenate((ref, np.zeros((1,), dtype=np.int8)))
        coord = np.array(g[-3::2], dtype=np.str).astype(np.float64)
        coord = np.concatenate((coord, np.zeros((1,), dtype=np.float64)))
        index_array[i-1] = ref
        coord_array[i-1] = coord
    else:
        atom_types[i] = g[:-6]
        index_array[i-1] = np.array(g[-6::2], dtype=np.str).astype(np.int8)
        coord_array[i-1] = np.array(g[-5::2], dtype=np.str).astype(np.float64)

zmat_re_pattern = grp_p(name_p)
for i in range(3):
    zmat_re_pattern += opnb_p(# optional non-binding
        wsr_p + grp_p(posint_p) + # ref in as a group
        wsr_p + grp_p(num_p) # ref value as a group
    )

zmat_res = re.compile(zmat_re_pattern)
def pull_zmat(txt,
             num_atoms = None,
             regex = zmat_res,
             num_header = 1
             ):
    """Pulls coordinates out of a zmatrix

    :param txt:
    :type txt:
    :param num_atoms:
    :type num_atoms:
    :param regex:
    :type regex:
    :return:
    :rtype:
    """

    if num_atoms is None:
        num_cur = 15 # probably more than we'll actually need...
        atom_types = [None]*num_cur
        index_array = np.zeros((num_cur-1, 3))
        coord_array = np.zeros((num_cur-1, 3))
        for i, match in enumerate(re.finditer(regex, txt)):
            if i == num_cur:
                atom_types.extend([None]*(2*num_cur))
                coord_array = np.concatenate(
                    (coord_array, np.zeros((2*num_cur, 6))),
                    axis=1
                )

            g = match.groups()
            if i == 0:
                g = g[:-6]
            elif i == 1:
                g = g[:-4]
            elif i == 2:
                g = g[:-2]
            # print(g)
            process_zzzz(i, g, atom_types, index_array, coord_array,  num_header=num_header)
        atom_types = atom_types[:i+1]
        index_array = index_array[:i]
        coord_array = coord_array[:i]

    else:
        atom_types = [None]*num_atoms
        index_array = np.zeros((num_atoms-1, 3))
        coord_array = np.zeros((num_atoms-1, 3))
        parse_iter = re.finditer(regex, txt)
        for i in range(num_atoms):
            match = next(parse_iter)
            g = match.groups()
            process_zzzz(i, g, atom_types, index_array, coord_array, num_header=num_header)

    return (atom_types, index_array, coord_array)


