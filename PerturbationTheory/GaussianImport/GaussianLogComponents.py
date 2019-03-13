"""This lists the types of readers and things available to the GaussianLogReader

"""

GuassianLogComponents = { } # we'll register on this bit by bit
# each registration should look like:

# GuassianLogComponents["Name"] = {
#     "TagStart" : start_tag, # starting delimeter for a block
#     "TagEnd"   : end_tag, # ending delimiter for a block None means apply the parser upon TagStart
#     "Parser"   : parser, # function that'll parse the returned list of blocks (for "List") or single block (for "Single)
#     "Mode"     : mode # "List" or "Single"
#
# }


########################################################################################################################
#
#                                           CartesianCoordinates
#

cartesian_coordinates_tags = (
""" ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------""",
""" ---------------------------------------------------------------------""",

)

cartesian_coordinates_parser = None


GuassianLogComponents["CartesianCoordinates"] = {
    "TagStart" : cartesian_coordinates_tags[0],
    "TagEnd"   : cartesian_coordinates_tags[1],
    "Parser"   : cartesian_coordinates_parser,
    "Mode"     : "List"
}