"""This lists the types of readers and things available to the GaussianLogReader

"""

from .ParserUtils import *

########################################################################################################################
#
#                                           GaussianLogComponents
#
#region GaussianLogComponents
GaussianLogComponents = { } # we'll register on this bit by bit
# each registration should look like:

# GaussianLogComponents["Name"] = {
#     "description" : string, # used for docmenting what we have
#     "tag_start"    : start_tag, # starting delimeter for a block
#     "tag_end"      : end_tag, # ending delimiter for a block None means apply the parser upon tag_start
#     "parser"      : parser, # function that'll parse the returned list of blocks (for "List") or single block (for "Single)
#     "mode"        : mode # "List" or "Single"
#
# }

########################################################################################################################
#
#                                           InputZMatrix
#

#region InputZMatrix
tag_start  = "Z-matrix:"
tag_end    = """     
 
"""
def parser(zmat):
    return zmat
mode       = "Single"

GaussianLogComponents["InputZMatrix"] = {
    "tag_start" : tag_start,
    "tag_end"   : tag_end,
    "parser"   : parser,
    "mode"     : mode
}

#endregion

########################################################################################################################
#
#                                           CartesianCoordinates
#

#region CartesianCoordinates

 # the region thing is just a PyCharm hack to collapse the boilerplate here... Could also have done 5000 files

cartesian_coordinates_tags = (
""" ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------""",
""" ---------------------------------------------------------------------""",

)

cartesian_re_c = re.compile(ws_p.join(["("+int_p+")"]*3)+ws_p+cart_p)
def cartesian_coordinates_parser(strs):
    num_sets = len(strs)
    strit = iter(strs)
    if num_sets>0:
        first = pull_xyz(next(strit), regex=cartesian_re_c)
        num_atoms = len(first[0])
        if num_sets>1:
            coord_array = np.zeros((num_sets, num_atoms, 3), dtype=np.float64)
            coord_array[0] = first[1]
            for i, xyz in enumerate(strs):
                if i>0:
                    coord_array[i] = pull_coords(xyz)
            coords = [first[0], coord_array]
        else:
            coords = [first[0], np.reshape(first[1], (1,) + first[1].shape)]
    else:
        coords = None

    return coords

GaussianLogComponents["CartesianCoordinates"] = {
    "tag_start" : cartesian_coordinates_tags[0],
    "tag_end"   : cartesian_coordinates_tags[1],
    "parser"   : cartesian_coordinates_parser,
    "mode"     : "List"
}

#endregion

########################################################################################################################
#
#                                           ZMatrices
#

#region ZMatrices
tag_start  = """Z-MATRIX (ANGSTROMS AND DEGREES)
   CD    Cent   Atom    N1       Length/X        N2       Alpha/Y        N3        Beta/Z          J
 ---------------------------------------------------------------------------------------------------"""
tag_end    = " ---------------------------------------------------------------------"

gaussian_zzz = op_p(posint_p) + opnb_p(wsr_p) + \
                op_p(posint_p) + opnb_p(wsr_p) + \
                grp_p(name_p+opnb_p(paren_p))
for i in range(3):
    gaussian_zzz += opnb_p( # optional non-binding
        wsr_p + grp_p(posint_p) + # ref int
        wsr_p + grp_p(num_p) + # z-mat value
        opnb_p(paren_p) # ignore parens that Gaussian puts at the end
    )
# print(gaussian_zzz)
gaussian_zzz_c = re.compile(gaussian_zzz)
def parser(strs):
    num_sets = len(strs)
    strit = iter(strs)
    if num_sets>0:
        zm = next(strit)
        # print(zm)
        first = pull_zmat(zm, regex=gaussian_zzz_c, num_header=3)
        # print(first)
        num_atoms = len(first[0])
        if num_sets>1:
            index_array = first[1]
            coord_array = np.zeros((num_sets, num_atoms-1, 3), dtype=np.float64)
            coord_array[0] = first[2]
            for i, zm in enumerate(strs):
                if i>0:
                    coord_array[i] = pull_zmat_coords(zm)
            coords = [first[0], index_array, coord_array]
        else:
            coords = [first[0], first[1], np.reshape(first[2], (1,) + first[2].shape)]
    else:
        coords = None

    return coords
mode       = "List"

GaussianLogComponents["ZMatrices"] = {
    "tag_start" : tag_start,
    "tag_end"   : tag_end,
    "parser"   : parser,
    "mode"     : mode
}

#endregion

########################################################################################################################
#
#                                           MullikenCharges
#

#region MullikenCharges
tag_start  = "Mulliken charges:"
tag_end    = "Sum of Mulliken charges"
def parser(charges):
    """Parses a Mulliken charges block"""
    return charges
mode       = "List"

GaussianLogComponents["MullikenCharges"] = {
    "tag_start" : tag_start,
    "tag_end"   : tag_end,
    "parser"   : parser,
    "mode"     : mode
}

#endregion

########################################################################################################################
#
#                                           MultipoleMoments
#

#region MultipoleMoments
tag_start  = "Dipole moment ("
tag_end    = " N-N="
def parser(moms):
    """Parses a multipole moments block"""
    return moms
mode       = "List"

GaussianLogComponents["MultipoleMoments"] = {
    "tag_start" : tag_start,
    "tag_end"   : tag_end,
    "parser"   : parser,
    "mode"     : mode
}

#endregion

########################################################################################################################
#
#                                           DipoleMoments
#

#region DipoleMoments
tag_start  = "Dipole moment ("
tag_end    = "Quadrupole moment ("
def parser(moms):
    """Parses a multipole moments block"""
    return moms
mode       = "List"

GaussianLogComponents["DipoleMoments"] = {
    "tag_start" : tag_start,
    "tag_end"   : tag_end,
    "parser"   : parser,
    "mode"     : mode
}

#endregion

#endregion

########################################################################################################################
#
#                                           GaussianLogDefaults
#
#region GaussianLogDefaults
GaussianLogDefaults = (
    "StartDateTime",
    "InputZMatrix",
    "ScanTable",
    "Blurb",
    "ComputerTimeElapsed",
    "EndDateTime"
)
#endregion

########################################################################################################################
#
#                                           GaussianLogOrdering
#
#region GaussianLogOrdering
# defines the ordering in a GaussianLog file
glk = ( # this must be sorted by what appears when
    "StartDateTime",
    "AtomPositions",
    "CartesianCoordinates",
    "CartesianCoordinateVectors",
    "MullikenCharges",
    "MultipoleMoments",
    "DipoleMoments",
    "QuadrupoleMoments",
    "OctapoleMoments",
    "HexadecapoleMoments",
    "HartreeFockEnergies",
    "MP2Energies",
    "InputZMatrix",
    "InputZMatrixVariables",
    "ZMatrices",
    "ZMatrixCoordinates",
    "ZMatrixCoordinateVectors",
    "ScanTable",
    "OptimizationScan",
    "Blurb",
    "ComputerTimeElapsed",
    "EndDateTime"
)
GaussianLogOrdering = { k:i for i, k in enumerate(glk) }
del glk
#endregion