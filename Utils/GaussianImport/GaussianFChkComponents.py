"""Defines components of an .fchk file that are already known and parseable"""

from .ParserUtils import *
from .FChkDerivatives import *

########################################################################################################################
#
#                                           FormattedCheckpointComponents
#
#region FormattedCheckpointComponents
FormattedCheckpointComponents = { } # we'll register on this bit by bit
# each registration should look like:

# FormattedCheckpointComponents["Name"] = parser

########################################################################################################################
#
#                                          Current cartesian coordinates
#

#region Current cartesian coordinates

def reformat(coords):
    ncoords = len(coords)
    return np.reshape(coords, (int(ncoords/3), 3))
FormattedCheckpointComponents["Current cartesian coordinates"] = reformat

#endregion

########################################################################################################################
#
#                                           Cartesian Force Constants
#

#region Cartesian Force Constants

FormattedCheckpointComponents["Cartesian Force Constants"] = FchkForceConstants

#endregion

########################################################################################################################
#
#                                           Cartesian 3rd/4th derivatives
#

#region Cartesian 3rd/4th derivatives

FormattedCheckpointComponents["Cartesian 3rd/4th derivatives"] = FchkForceDerivatives

#endregion

########################################################################################################################
#
#                                           Vib-Modes
#

#region Vib-Modes

def split_vib_modes(mcoeffs):
    """Pulls the mode vectors from the coeffs
    There should be 3N-6 modes where each vector is 3N long so N = (1 + sqrt(1 + l/9))

    :param mcoeffs:
    :type mcoeffs:
    :return:
    :rtype:
    """
    l = len(mcoeffs)
    n = int(1 + np.sqrt(1 + l/9))
    return np.reshape(mcoeffs, (3*n-6, 3*n))
FormattedCheckpointComponents["Vib-Modes"] = split_vib_modes

#endregion

#endregion


########################################################################################################################
#
#                                           CommonNames
#

#region CommonNames

FormattedCheckpointCommonNames = {

    "Atomic numbers": "AtomicNumbers",
    "Current cartesian coordinates":"Coordinates",
    "Cartesian Force Constants" : "ForceConstants",
    "Cartesian 3rd/4th derivatives" : "ForceDerivatives",
    "Vib-E2" : "VibrationalEnergies",
    "Vib-Modes" : "VibrationalModes",
    "Vib-AtMass" : "VibrationalMasses",
    "Real atomic weights" : "AtomicMasses"

}

#endregion