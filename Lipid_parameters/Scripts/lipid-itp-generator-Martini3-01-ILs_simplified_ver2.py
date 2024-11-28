#!/usr/bin/env python3
"""lipid-itp-generator-Martini3 creates a customized Martini lipid topologies for Martini 3 ionizable lipids, use lipid-itp-generator-Martini3.py -h for description"""

__author__  = "Helgi I. Ingolfsson and Kasper Busk Pedersen, Tsjerk A. Wassenaar, and Lisbeth Kjølbye"
__status__  = "Development"
__version__ = "M3.l01 Ionizable lipids"
__email__   = ""

import sys,math

# Very simple option class
class Option:
    def __init__(self,func=str,num=1,default=None,description=""):
        self.func        = func
        self.num         = num
        self.value       = default
        self.description = description
    def __nonzero__(self):
        if self.func == bool:
            return self.value != False
        return bool(self.value)
    def __str__(self):
        return self.value and str(self.value) or ""
    def setvalue(self,v):
        if len(v) == 1:
            self.value = self.func(v[0])
        else:
            self.value = [ self.func(i) for i in v ]


# Description
desc = """
This script creates a customized Martini 3 lipid topology based on the head and
tail specification strings provided. The topology for tails follows the standard Martini 3.0 lipid
definitions as described in:
  P.C.T. Souza et al. Martini 3: a general purpose force field for coarse-grained molecular dynamics,
  Nat. Methods; 2021. doi: 10.1038/s41592-021-01098-3
This scripts is based on the Martini 2 lipid topology builder script lipid-martini-itp-v06.py
and described in:
  T.A. Wassenaar, H.I. Ingolfsson, R.A. Bockmann, D.P. Tieleman, S.J. Marrink. Computational
  lipidomics with insane: a versatile tool for generating custom membranes for molecular simulations.
  JCTC, 150410125128004, 2015. doi:10.1021/acs.jctc.5b00209


WARNING:
  - This script can generate topologies for numerous lipids many of which are unrealistic
  and untested, please use with discretion.


The lipid descriptions supported are as follows:

Heads (-alhead):
  Please provide a list of lipid head beads. The left most bead will be on top, they are
  connected in a sequence from left to right and the right most bead is connected to the
  first bead in the linker. NB the first tail bead can only be C1, a saturated bead. 

  head bead types supported:
@TODO    SD = Tertiary amine with short link to alcohol         - bead SN3a,  charge 0
@TODO    SP = Tertiary amine with short link to alcohol         - bead SQ2p,  charge +1
@TODO    AD = Tertiary amine with long link to alcohol          - bead SN3a,  charge 0
@TODO    AP = Tertiary amine with long link to alcohol          - bead SQ2p,  charge +1

    M3 = Tertiary amine with one ester as linker                         - bead SN3a,  charge 0
    MP = Tertiary amine with one ester as linker                         - bead SQ2p,  charge +1
    K2 = Tertiary amine linked to one ketal                              - bead SN3a,  charge 0
    KP = Tertiary amine linked to one ketal                              - bead SQ2p,  charge +1
    B1 = Tertiary amine with one ketal as linker                         - bead SN3a,  charge 0
    BP = Tertiary amine with one ketal as linker                         - bead SQ2p,  charge +1
    DT = Tertiary amine with ether as linkers                           - bead SN3a,  charge 0
    DP = Tertiary amine with ether as linkers                           - bead SQ2p,  charge +1
    NE = Tertiary amine with ester as linkers                           - bead SN3a,  charge 0
    NP = Tertiary amine with ester as linkers                           - bead SQ2p,  charge +1
    AE = Tertiary amine with carbon link to ester branching              - bead SN3a,  charge 0
    AP = Tertiary amine with carbon link to ester branching              - bead SN3a,  charge 0
    PI = Piperazine with carbon link to ester branching                  - bead N6a,  charge 0
    PP = Piperazine with carbon link to ester branching                  - bead Q2p,  charge +2
    S1 = Alcohol link to tertiary amide link to ester                    - bead SN3a,  charge 0
    SP = Alcohol link to tertiary amide link to ester                    - bead SQ2p,  charge +1
    Z1 = Larger Alcohol link to tertiary amide link to ester             - bead SN3a,  charge 0
    ZP = Larger Alcohol link to tertiary amide link to ester             - bead SQ2p,  charge +1


Tails (-altail):
  One lipid tail definition should be provided for each linker, separated with a space;
  extra spaces are ignored. Each tail can have an arbitrary number of tail beads.

  tail bead types supported:
    C = corresponding roughly to a linear combination of x4 CH2 groups (CH2-CH2-CH2-CH2).
        Represented with a C1 bead.
    c = corresponding roughly to a linear combination of x2 CH2 groups (CH2-CH2).
        Represented with a SC1 bead.
    D = corresponding roughly to a linear combination of x2 CH2 and x2 CH groups
        (CH2-CH=CH-CH2), where the double bound is in a cis bond. Represented with a
        C4h bead.
    F = corresponding roughly to a linear combination carbon groups with more than
        one double bond (noramlly 1.5). Represented with a C5h bead.
@TODO check T.
    T = Same as D except with a trans double bond. The angle [X / T / X] is set with
        equilibrium bond angle Theta_a = 180 and K_a = 45. Represented with a C3 bead.

@TODO update
  Examples of tails (first bead of each tail has to be C):
    Lyso tails:
    "- CCCC       " - C16-18:0 Lyso
    "- CCCCC      " - C20-22:0 Lyso
    "- CCCCCC     " - C24-26:0 Lyso
    Saturated tails:
    "CC     CC    " - C08-10:0 - diOctanoyl or diDecanoyl
    "CCC    CCC   " - C12-14:0 - diLauric acid or diMyristoyl
    "CCCC   CCCC  " - C16-18:0 - diPalmitic acid or diStearoyl
    "CCCCC  CCCCC " - C20-22:0 - diArachidoyl or diBehenoyl
    "CCCCCC CCCCCC" - C24-26:0 - diLignoceroyl or diHexacosanoyl
    Unsaturated tails:
    "CDC    CDC   " - C14:1(9c) - diMyristoleoyl
    "CDCC   CDCC  " - C16:1(9c) - diPalmitoleoyl / C18:1(9c) - diOleoyl
    "CDDC   CDDC  " - C18:2(9c,12c) - diLinoleoyl
    "CCDCC  CCDCC " - C20:1(11c) - diGondic acid / C22:1(11c) - diErucoyl
    "CCCDCC CCCDCC" - C24:1(15c) - diNervonoyl
    Mixed tails:
    "CCCC   CCC   " - C14:0/C16:0 - MP / C14:0/C18:0 - MS
    "CDCC   CCCC  " - C16:0/C18:1(9c) - PO / C18:0/C18:1(9c) - SO
    "CDDC   CCCC  " - C16:0/C18:2(9c,12c) / C18:0/C18:2(9c,12c)

    NOTE: the first tail (tail A) is connected to linker 1 closer to head (this is sn-2 for GLY linker lipids), which is reverse order
    compared to how regular lipid names are written. The second tail is tail B (for GLY linker lipids this is sn-1)

@Update
Use:
  ./lipid-itp-generator-Martini3-01.py -alhead 'C P' -allink 'G G' -altail "CDCC cCCC" -alname POPC -o POPC-lipid.itp
"""

# Options
options = [
"""
Options:""",
("-o",       Option(str,    1,        "Martini-lipid.itp", "Output speciffic Martini lipid topology")),
("-alname",  Option(str,    1,        "POPC", "Four letter lipid name")),
("-alhead",  Option(str,    1,        "C P", "Lipid heads, see description")),
("-altail",  Option(str,    1,        "CDCC cCCC", "Lipid tails, see description")),
("-name",    Option(str,    1,        "POPC", "A common name of the lipid, only use in comments")),
("-desc",    Option(str,    1,        "This is a ...", "A general description of what the FF is / represents, only use in comments")),
("-keyw",    Option(str,    1,        "", "List of keywords, only use in comments")),
("-parm",    Option(str,    1,        "Was modeled on ...", "Fow the FF was parameterized, only use in comments")),
("-refs",    Option(str,    1,        "", "List of references for the FF, only use in comments")),
("-crea",    Option(str,    1,        "", "FF created on, only use in comments")),
("-auth",    Option(str,    1,        "", "FF author, only use in comments")),
("-modi",    Option(str,    1,        "", "List of modifications to the FF, only use in comments")),
("-area",    Option(str,    1,        "", "Reference area per lipid, only use in comments")),
("-warn",    Option(str,    1,        "", "Warning(s)/Note(s) for the FF, only use in comments"))
          ]

# Define supported lipid head beads
# Lists all supported head bead types. One letter name mapped to type, atom name and charge
headMapp = {
    # beadtype, beadname, charge, ffbonded name
    "C":  ['Q1',   'NC3', '1.0',  'NC3'],  # NC3 = Choline
    "E":  ['Q4p',  'NH3', '1.0',  'NH3'],  # NH3 = Ethanolamine
    "G":  ['P4',   'GL0', '0.0',  'GL0'],  # GL0 = Glycerol (used for PG lipid head group)
    "S":  ['P5',   'CNO', '0.0',  'CNO'],  # CNO = Serine  - x1 bead PS WARNING not been updated here - used x2 bead version (PS1 PS2)
  "PS1":  ['SP2q', 'PS1', '-0.6', 'PS1'],  # PS1 is bead one of two bead PS represents a COO group
  "PS2":  ['Q4p',  'PS2', '0.6',  'PS2'],  # PS2 is bead two of two bead PS represents a NH3 group
    "P":  ['Q5',   'PO4', '-1.0', 'PO4'],  # PO4 = Phosphate
    "O":  ['Q5',   'PO4', '-2.0', 'PO4'],  # PO4 = Phosphate (one bond x2 charges can be used e.g. when making unprotonated PA lipids)
  "COH":  ['TN6',  'COH', '0.0',  'COH'],  # COH = Cappding bead for top of diacylglycerols and ceramides
    }


tailMapp = {
    # beadtype, beadname, charge, ffbonded bondname and angle/dih name
    "C":  ['C1',   'C??', '0',  'C1',  'C1'],  # C = straight chain
    "c":  ['SC1',  'C??', '0', 'SC1',  'C1'],  # c = short straight chain
    "D":  ['C4h',  'D??', '0',  'C4',  'C4'],  # D = chain with double bond
    "F":  ['C5h',  'D??', '0',  'C4',  'C4'],  # F = chain with more than one double bond (normaly 1.5)
    "T":  ['C4h',  'T??', '0',  'C4',  'C1'],  # T = chain with trans double bond - normal sized bead - only used in shingosine top bead (SM and CER)
    "t":  ['SC4h', 'T??', '0', 'SC4',  'C1'],  # t = chain with trans double bond - small sized bead - only used in shingosine top bead (SM and CER)
    }

# Get arguments
args = sys.argv[1:]

# Print help
if '-h' in args or '--help' in args:
    print("\n", __file__)
    print(desc)
    for thing in options:
        print(type(thing) != str and "%10s  %s"%(thing[0],thing[1].description) or thing)
    print()
    sys.exit()

# Convert the option list to a dictionary, discarding all comments
options = dict([i for i in options if not type(i) == str])

# Process the command line
while args:
    ar = args.pop(0)
    options[ar].setvalue([args.pop(0) for i in range(options[ar].num)])

# Get ouput .itp file name
itpFileName  = options["-o"].value

# Get lipid description
lipidHead  = options["-alhead"].value
lipidTail  = options["-altail"].value
lipidName = options["-alname"].value

lipidCommonName = options["-name"].value
lipidDesc = options["-desc"].value
lipidParm = options["-parm"].value
if lipidCommonName==None or lipidDesc==None or lipidParm==None:
    print("You have to provide a common name, description and list how the FF was parameterized.", file=sys.stderr)
    sys.exit()
lCharge = 0  # Update when adding charged beads

progString = "The Martini lipid itp generator version " + __version__ + "  Args are: -alname %s -alhead '%s' -altail '%s'" % (lipidName, lipidHead, lipidTail)
print(f"{progString} to file {itpFileName}")

headsArray = lipidHead.split()
tailsArray = lipidTail.split()
#print (tailsArray)
#print (tailsArray[0][:2])

#drawingArray = [[]]
#drawingHeadsize = 0
bondsArray = []
anglesArray = []
beadArray = []
dihedralsArray = []
constraintsArray = []
exclusionsArray = []

# If speciall head insert now all beads, bonds, angles, dihedrals, constreints etc
index = 1


### ILS ####


if len(headsArray)>0 and headsArray[0]=='Z1':
    beadArray.append([1,  'P1  ',     1,   lipidName,    'OH',         1,      0,     ''])
    beadArray.append([2,  'SN3a',     1,   lipidName,    'N1',         2,      0,     ''])
    beadArray.append([3,  'C2  ',     1,   lipidName,    'CA',         3,      0,     ''])
    beadArray.append([4,  'SN4a',     1,   lipidName,    'GLA',        4,      0,     ''])
    beadArray.append([5,  'C2  ',     1,   lipidName,    'CX ',        5,      0,     '']) 
    index += 5
    lCharge += 0.0 # Keep track of overall lipid charge
    bondsArray.append([-1, 'Headgroup and linker bonds'])
    bondsArray.append([1,      2,      1,  '0.52',   '2500',  'OH-N1'])
    bondsArray.append([2,      3,      1,  '0.52',   '2500',  'N1-CA'])
    bondsArray.append([3,      4,      1,  '0.37',   '5000',  'CA-GLA'])
    bondsArray.append([4,      5,      1,  '0.37',   '5000',  'GLA-CX'])
    bondsArray.append([5 ,     6,      1,  '0.47',   '5000',  'CX-C1A'])
    bondsArray.append([5 ,   len(tailsArray[0])+index  , 1,  '0.47','5000', 'CX-C1B'])
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])

    anglesArray.append([-1, 'Orient the headgroup and linker'])
    anglesArray.append([1,      2,      3,      2,   '115',   '32',  'OH-N1-CA'])
    anglesArray.append([4,   5,   len(tailsArray[0])+index,   2,  '100', '25', 'GLA-CX-C1B'])
    anglesArray.append([4,   5,   6,   2,  '100', '25', 'GLA-CX-C1A'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])
    if tailsArray[0][:2] == 'CC':
        anglesArray.append([5,       index,       1+index,   2,     '155',      '25', 'CX-C1A-C2A']) 
    elif tailsArray[0][:2] == 'CD':
        anglesArray.append([5,       index,       1+index,   2,     '140',      '15', 'CX-C1A-D2A'])
    elif tailsArray[0][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[0][:2] == 'DD':
        print ('Not possible yet')
    if tailsArray[1][:2] == 'CC':
        anglesArray.append([5,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '155',      '25', 'CX-C1B-C2B']) 
    elif tailsArray[1][:2] == 'CD':
        anglesArray.append([5,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '140',      '15', 'CX-C1B-D2B'])
    elif tailsArray[1][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[1][:2] == 'DD':
        print ('Not possible yet')
    exclusionsArray.append([-2, '1 2 3'])
    exclusionsArray.append([-2, '6 {}'.format((len(tailsArray[0])+index))])
    exclusionsArray.append([-2, '[ pairs ]'])
    exclusionsArray.append([-2, '6 {} 1 4.100000e-01    2.350000e+00'.format(len(tailsArray[0])+index)])





if len(headsArray)>0 and headsArray[0]=='ZP':
    beadArray.append([1,  'P1  ',     1,   lipidName,    'OH',         1,      0,     ''])
    beadArray.append([2,  'SQ2p',     1,   lipidName,    'NP',         2,      1,     ''])
    beadArray.append([3,  'C2  ',     1,   lipidName,    'CA',         3,      0,     ''])
    beadArray.append([4,  'SN4a',     1,   lipidName,    'GLA',        4,      0,     ''])
    beadArray.append([5,  'C2  ',     1,   lipidName,    'CX ',        5,      0,     '']) 
    index += 5
    lCharge += 1.0 # Keep track of overall lipid charge
    bondsArray.append([-1, 'Headgroup and linker bonds'])
    bondsArray.append([1,      2,      1,  '0.52',   '2500',  'OH-NP'])
    bondsArray.append([2,      3,      1,  '0.52',   '2500',  'NP-CA'])
    bondsArray.append([3,      4,      1,  '0.37',   '5000',  'CA-GLA'])
    bondsArray.append([4,      5,      1,  '0.37',   '5000',  'GLA-CX'])
    bondsArray.append([5 ,     6,      1,  '0.47',   '5000',  'CX-C1A'])
    bondsArray.append([5 ,   len(tailsArray[0])+index  , 1,  '0.47','5000', 'CX-C1B'])
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])

    anglesArray.append([-1, 'Orient the headgroup and linker'])
    anglesArray.append([1,      2,      3,      2,   '115',   '30',  'OH-N1-CA'])
    anglesArray.append([4,   5,   len(tailsArray[0])+index,   2,  '100', '25', 'GLA-CX-C1B'])
    anglesArray.append([4,   5,   6,   2,  '100', '25', 'GLA-CX-C1A'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])
    if tailsArray[0][:2] == 'CC':
        anglesArray.append([5,       index,       1+index,   2,     '155',      '25', 'CX-C1A-C2A']) 
    elif tailsArray[0][:2] == 'CD':
        anglesArray.append([5,       index,       1+index,   2,     '140',      '15', 'CX-C1A-D2A'])
    elif tailsArray[0][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[0][:2] == 'DD':
        print ('Not possible yet')
    if tailsArray[1][:2] == 'CC':
        anglesArray.append([5,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '155',      '25', 'CX-C1B-C2B']) 
    elif tailsArray[1][:2] == 'CD':
        anglesArray.append([5,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '140',      '15', 'CX-C1B-D2B'])
    elif tailsArray[1][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[1][:2] == 'DD':
        print ('Not possible yet')
    exclusionsArray.append([-2, '1 2 3'])
    exclusionsArray.append([-2, '6 {}'.format((len(tailsArray[0])+index))])
    exclusionsArray.append([-2, '[ pairs ]'])
    exclusionsArray.append([-2, '6 {} 1 4.100000e-01    2.350000e+00'.format(len(tailsArray[0])+index)])





if len(headsArray)>0 and headsArray[0]=='S1':
    beadArray.append([1,  'TP1 ',      1    , lipidName,    'OH ',      1,      0,     ''])
    beadArray.append([2,  'SN4a',      1    , lipidName,    'N1 ',      2,      0,     ''])
    beadArray.append([3,  'SC2 ',      1    , lipidName,    'CA ',      3,      0,     ''])
    beadArray.append([4,  'SN4a',      1    , lipidName,    'GLA',      4,      0,     ''])
    beadArray.append([5,  'C2  ',      1    , lipidName,    'CX ',      5,      0,     '']) 
    index += 5
    lCharge += 0.0 # Keep track of overall lipid charge
    bondsArray.append([-1, 'Headgroup and linker bonds'])
    bondsArray.append([1,      2,      1,  '0.34',   '5000',  'OH-N1'])
    bondsArray.append([2,      3,      1,  '0.50',   '5000',  'N1-CA'])
    bondsArray.append([3,      4,      1,  '0.37',   '5000',  'CA-GLA'])
    bondsArray.append([4,      5,      1,  '0.37',   '5000',  'GLA-CX'])
    bondsArray.append([5 ,     6,      1,  '0.47',   '5000',  'CX-C1B'])
    bondsArray.append([5 ,   len(tailsArray[0])+index  , 1,  '0.47','5000', 'CX-C1B'])
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])

    anglesArray.append([-1, 'Orient the headgroup and linker'])
    anglesArray.append([1,      2,      3,      2,   '118',   '15',  'OH-N1-CA'])
    anglesArray.append([4,   5,   len(tailsArray[0])+index,   2,  '100', '25', 'GLA-CX-C1B'])
    anglesArray.append([4,   5,   6,   2,  '100', '25', 'GLA-CX-C1A'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])
    if tailsArray[0][:2] == 'CC':
        anglesArray.append([5,       index,       1+index,   2,     '155',      '25', 'CX-C1A-C2A']) 
    elif tailsArray[0][:2] == 'CD':
        anglesArray.append([5,       index,       1+index,   2,     '140',      '15', 'CX-C1A-D2A'])
    elif tailsArray[0][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[0][:2] == 'DD':
        print ('Not possible yet')
    if tailsArray[1][:2] == 'CC':
        anglesArray.append([5,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '155',      '25', 'CX-C1B-C2B']) 
    elif tailsArray[1][:2] == 'CD':
        anglesArray.append([5,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '140',      '15', 'CX-C1B-D2B'])
    elif tailsArray[1][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[1][:2] == 'DD':
        print ('Not possible yet')
    exclusionsArray.append([-2, '6 {}'.format((len(tailsArray[0])+index))])
    exclusionsArray.append([-2, '[ pairs ]'])
    exclusionsArray.append([-2, '6 {} 1 4.100000e-01    2.350000e+00'.format(len(tailsArray[0])+index)])




if len(headsArray)>0 and headsArray[0]=='SP':
    beadArray.append([1,  'TP1 ',      1    , lipidName,    'OH ',      1,      0,     ''])
    beadArray.append([2,  'SQ2p',      1    , lipidName,    'NP ',      2,      1,     ''])
    beadArray.append([3,  'SC2 ',      1    , lipidName,    'CA ',      3,      0,     ''])
    beadArray.append([4,  'SN4a',      1    , lipidName,    'GLA',      4,      0,     ''])
    beadArray.append([5,  'C2  ',      1    , lipidName,    'CX ',      5,      0,     '']) 
    index += 5
    lCharge += 1.0 # Keep track of overall lipid charge
    bondsArray.append([-1, 'Headgroup and linker bonds'])
    bondsArray.append([1,      2,      1,  '0.34',   '5000',  'OH-NP'])
    bondsArray.append([2,      3,      1,  '0.47',   '5000',  'NP-CA'])
    bondsArray.append([3,      4,      1,  '0.37',   '5000',  'CA-GLA'])
    bondsArray.append([4,      5,      1,  '0.37',   '5000',  'GLA-CX'])
    bondsArray.append([5 ,     6,      1,  '0.47',   '5000',  'CX-C1B'])
    bondsArray.append([5 ,   len(tailsArray[0])+index  , 1,  '0.47','5000', 'CX-C1B'])
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])

    anglesArray.append([-1, 'Orient the headgroup and linker'])
    anglesArray.append([1,      2,      3,      2,   '118',   '15',  'OH-N1-CA'])
    anglesArray.append([4,   5,   len(tailsArray[0])+index,   2,  '100', '25', 'GLA-CX-C1B'])
    anglesArray.append([4,   5,   6,   2,  '100', '25', 'GLA-CX-C1A'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])
    if tailsArray[0][:2] == 'CC':
        anglesArray.append([5,       index,       1+index,   2,     '155',      '25', 'CX-C1A-C2A']) 
    elif tailsArray[0][:2] == 'CD':
        anglesArray.append([5,       index,       1+index,   2,     '140',      '15', 'CX-C1A-D2A'])
    elif tailsArray[0][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[0][:2] == 'DD':
        print ('Not possible yet')
    if tailsArray[1][:2] == 'CC':
        anglesArray.append([5,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '155',      '25', 'CX-C1B-C2B']) 
    elif tailsArray[1][:2] == 'CD':
        anglesArray.append([5,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '140',      '15', 'CX-C1B-D2B'])
    elif tailsArray[1][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[1][:2] == 'DD':
        print ('Not possible yet')
    exclusionsArray.append([-2, '6 {}'.format((len(tailsArray[0])+index))])
    exclusionsArray.append([-2, '[ pairs ]'])
    exclusionsArray.append([-2, '6 {} 1 4.100000e-01    2.350000e+00'.format(len(tailsArray[0])+index)])


if len(headsArray)>0 and headsArray[0]=='PI':
    beadArray.append([1,  'N4a',       1    , lipidName,    'N1 ',      1,      0 ,    '']) 
    beadArray.append([2,   'N4a  ',     1    , lipidName,    'N2 ',      2,      0 ,    '']) 
    beadArray.append([3,   'C2  ',      1    , lipidName,    'CN' ,      3,      0 ,    '']) 
    beadArray.append([4,   'SN4a',      1    , lipidName,    'GLA',      4,      0 ,    '']) 
    beadArray.append([5,   'C2  ',      1    , lipidName,    'CX' ,      5,      0 ,    '']) 
    index += 5
    lCharge += 0.0 # Keep track of overall lipid charge
    bondsArray.append([-1, 'Headgroup and linker bonds'])
    bondsArray.append([1,      2,      1,  '0.37',   '2500',  'N1-N2'])
    bondsArray.append([2,      3,      1,  '0.45',   '2500',  'N2-CN'])
    bondsArray.append([3,      4,      1,  '0.45',   '2500',  'CN-GLA'])
    bondsArray.append([4,      5,      1,  '0.37',   '5000',  'GLA-CX'])
    bondsArray.append([5 ,     6,      1,  '0.47',   '5000',  'CX-C1A'])
    bondsArray.append([5 ,   len(tailsArray[0])+index  , 1,  '0.47','5000', 'CX-C1B'])
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])

    anglesArray.append([-1, 'Orient the headgroup and linker'])
    anglesArray.append([1,      2,      3,      2,   '150',   '35',  'N1-N2-CN'])
    anglesArray.append([3,      4,      5,      2,   '120',   '35',  'CN-GLA-CX'])
    anglesArray.append([4,   5,   len(tailsArray[0])+index,   2,  '100', '25', 'GLA-CX-C1B'])
    anglesArray.append([4,   5,   6,   2,  '100', '25', 'GLA-CX-C1A'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])
    if tailsArray[0][:2] == 'CC':
        anglesArray.append([5,       index,       1+index,   2,     '155',      '25', 'CX-C1A-C2A']) 
    elif tailsArray[0][:2] == 'CD':
        anglesArray.append([5,       index,       1+index,   2,     '140',      '15', 'CX-C1A-D2A'])
    elif tailsArray[0][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[0][:2] == 'DD':
        print ('Not possible yet')
    if tailsArray[1][:2] == 'CC':
        anglesArray.append([5,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '155',      '25', 'CX-C1B-C2B']) 
    elif tailsArray[1][:2] == 'CD':
        anglesArray.append([5,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '140',      '15', 'CX-C1B-D2B'])
    elif tailsArray[1][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[1][:2] == 'DD':
        print ('Not possible yet')
    exclusionsArray.append([-2, '6 {}'.format(len(tailsArray[0])+index)])
    exclusionsArray.append([-2, '[ pairs ]'])
    exclusionsArray.append([-2, '6 {} 1 4.100000e-01    2.350000e+00'.format(len(tailsArray[0])+index)])




if len(headsArray)>0 and headsArray[0]=='PP':
    beadArray.append([1,  'Q2p',       1    , lipidName,    'NP1 ',      1,      1 ,    '']) 
    beadArray.append([2,  'Q2p  ',     1    , lipidName,    'NP2 ',      2,      1 ,    '']) 
    beadArray.append([3,  'C3  ',      1    , lipidName,    'CN' ,      3,      0 ,    '']) 
    beadArray.append([4,  'SN4a',      1    , lipidName,    'GLA',      4,      0 ,    '']) 
    beadArray.append([5,  'C2  ',      1    , lipidName,    'CX' ,      5,      0 ,    '']) 
    index += 5
    lCharge += 2.0 # Keep track of overall lipid charge
    bondsArray.append([-1, 'Headgroup and linker bonds'])
    bondsArray.append([1,      2,      1,  '0.37',   '2500',  'NP1-NP2'])
    bondsArray.append([2,      3,      1,  '0.45',   '2500',  'NP2-CN'])
    bondsArray.append([3,      4,      1,  '0.45',   '2500',  'CN-GLA'])
    bondsArray.append([4,      5,      1,  '0.37',   '5000',  'GLA-CX'])
    bondsArray.append([5 ,     6,      1,  '0.47',   '5000',  'CX-C1A'])
    bondsArray.append([5 ,   len(tailsArray[0])+index  , 1,  '0.47','5000', 'CX-C1B'])
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])

    anglesArray.append([-1, 'Orient the headgroup and linker'])
    anglesArray.append([1,      2,      3,      2,   '150',   '35',  'N1-N2-CN'])
    anglesArray.append([3,      4,      5,      2,   '120',   '35',  'CN-GLA-CX'])
    anglesArray.append([4,   5,   len(tailsArray[0])+index,   2,  '100', '25', 'GLA-CX-C1B'])
    anglesArray.append([4,   5,   6,   2,  '100', '25', 'GLA-CX-C1A'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])
    if tailsArray[0][:2] == 'CC':
        anglesArray.append([5,       index,       1+index,   2,     '155',      '25', 'CX-C1A-C2A']) 
    elif tailsArray[0][:2] == 'CD':
        anglesArray.append([5,       index,       1+index,   2,     '140',      '15', 'CX-C1A-D2A'])
    elif tailsArray[0][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[0][:2] == 'DD':
        print ('Not possible yet')
    if tailsArray[1][:2] == 'CC':
        anglesArray.append([5,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '155',      '25', 'CX-C1B-C2B']) 
    elif tailsArray[1][:2] == 'CD':
        anglesArray.append([5,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '140',      '15', 'CX-C1B-D2B'])
    elif tailsArray[1][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[1][:2] == 'DD':
        print ('Not possible yet')
    exclusionsArray.append([-2, '6 {}'.format(len(tailsArray[0])+index)])
    exclusionsArray.append([-2, '[ pairs ]'])
    exclusionsArray.append([-2, '6 {} 1 4.100000e-01    2.350000e+00'.format(len(tailsArray[0])+index)])


if len(headsArray)>0 and headsArray[0]=='AE':
    beadArray.append([1,    'SN3a',       1    , lipidName,    'N1 ',      1,      0 ,    '']) 
    beadArray.append([2,    'C2  ',       1    , lipidName,    'CN ',      2,      0 ,    '']) 
    beadArray.append([3,    'SN4a',       1    , lipidName,    'GLA',      3,      0 ,    '']) 
    beadArray.append([4,    'C2  ',       1    , lipidName,    'CX' ,      4,      0 ,    '']) 
    index += 4
    lCharge += 0.0 # Keep track of overall lipid charge
    bondsArray.append([-1, 'Headgroup and linker bonds'])
    bondsArray.append([1,      2,      1,  '0.47',   '2500',  'N1-CN'])
    bondsArray.append([2,      3,      1,  '0.42',   '2500',  'CN-GLA'])
    bondsArray.append([3,      4,      1,  '0.37',   '5000',  'GLA-CX'])
    bondsArray.append([4 ,     5,      1,  '0.47',   '5000',  'CX-C1A'])
    bondsArray.append([4 ,   len(tailsArray[0])+index  , 1,  '0.47','5000', 'CX-C1B'])
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])

    anglesArray.append([-1, 'Orient the headgroup and linker'])
    anglesArray.append([1,      2,      3,      2,   '140',   '15',  'N1-CN-GLA'])
    anglesArray.append([3,      4,      5,      2,   '120',   '35',  'CN-GLA-CX'])
    anglesArray.append([3,      4,      6,      2,   '100',   '25',  'GLA-CX-C1A'])
    anglesArray.append([3,      4,   len(tailsArray[0])+index,   2,  '100', '25', 'GLA-CX-C1B'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])
    if tailsArray[0][:2] == 'CC':
        anglesArray.append([4,       index,       1+index,   2,     '155',      '25', 'CX-C1A-C2A']) 
    elif tailsArray[0][:2] == 'CD':
        anglesArray.append([4,       index,       1+index,   2,     '140',      '15', 'CX-C1A-D2A'])
    elif tailsArray[0][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[0][:2] == 'DD':
        print ('Not possible yet')
    if tailsArray[1][:2] == 'CC':
        anglesArray.append([4,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '155',      '25', 'CX-C1B-C2B']) 
    elif tailsArray[1][:2] == 'CD':
        anglesArray.append([4,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '140',      '15', 'CX-C1B-D2B'])
    elif tailsArray[1][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[1][:2] == 'DD':
        print ('Not possible yet')
    exclusionsArray.append([-2, '5 {}'.format(len(tailsArray[0])+index)])
    exclusionsArray.append([-2, '[ pairs ]'])
    exclusionsArray.append([-2, '5 {} 1 4.100000e-01    2.350000e+00'.format(len(tailsArray[0])+index)])



if len(headsArray)>0 and headsArray[0]=='AP':
    beadArray.append([1,  'SQ2p',       1    , lipidName,    'NP ',      1,      1 ,    '']) 
    beadArray.append([2,  'C2  ',       1    , lipidName,    'CN ',      2,      0 ,    '']) 
    beadArray.append([3,  'SN4a',       1    , lipidName,    'GLA',      3,      0 ,    '']) 
    beadArray.append([4,  'C2  ',       1    , lipidName,    'CX' ,      4,      0 ,    '']) 
    index += 4
    lCharge += 1.0 # Keep track of overall lipid charge
    bondsArray.append([-1, 'Headgroup and linker bonds'])
    bondsArray.append([1,      2,      1,  '0.47',   '2500',  'NP-CN'])
    bondsArray.append([2,      3,      1,  '0.42',   '2500',  'CN-GLA'])
    bondsArray.append([3,      4,      1,  '0.37',   '5000',  'GLA-CX'])
    bondsArray.append([4 ,     5,      1,  '0.47',   '5000',  'CX-C1A'])
    bondsArray.append([4 ,   len(tailsArray[0])+index  , 1,  '0.47','5000', 'CX-C1B'])
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])

    anglesArray.append([-1, 'Orient the headgroup and linker'])
    anglesArray.append([1,      2,      3,      2,   '140',   '15',  'NP-CN-GLA'])
    anglesArray.append([3,      4,      5,      2,   '120',   '35',  'CN-GLA-CX'])
    anglesArray.append([3,      4,      6,      2,   '100',   '25', 'GLA-CX-C1A'])
    anglesArray.append([3,      4,   len(tailsArray[0])+index,   2,  '100', '25', 'GLA-CX-C1B'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])
    if tailsArray[0][:2] == 'CC':
        anglesArray.append([4,       index,       1+index,   2,     '155',      '25', 'CX-C1A-C2A']) 
    elif tailsArray[0][:2] == 'CD':
        anglesArray.append([4,       index,       1+index,   2,     '140',      '15', 'CX-C1A-D2A'])
    elif tailsArray[0][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[0][:2] == 'DD':
        print ('Not possible yet')
    if tailsArray[1][:2] == 'CC':
        anglesArray.append([4,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '155',      '25', 'CX-C1B-C2B']) 
    elif tailsArray[1][:2] == 'CD':
        anglesArray.append([4,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '140',      '15', 'CX-C1B-D2B'])
    elif tailsArray[1][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[1][:2] == 'DD':
        print ('Not possible yet')
    exclusionsArray.append([-2, '5 {}'.format(len(tailsArray[0])+index)])
    exclusionsArray.append([-2, '[ pairs ]'])
    exclusionsArray.append([-2, '5 {} 1 4.100000e-01    2.350000e+00'.format(len(tailsArray[0])+index)])



if len(headsArray)>0 and headsArray[0]=='M3':
    beadArray.append([1, 'TN3a', 1, lipidName,  'N1' , 1, 0,    ''])
    beadArray.append([2, 'TC2' , 1, lipidName,  'CN' , 2, 0,    ''])
    beadArray.append([3, 'SN4a', 1, lipidName,  'GLA', 3, 0,    ''])
    beadArray.append([4, 'SC2' , 1, lipidName,  'CX' , 4, 0, 'Branching bead to tail'])
    index += 4
    lCharge += 0.0 # Keep track of overall lipid charge

    bondsArray.append([-1, 'Headgroup and linker bonds'])

    bondsArray.append([1 ,   2  , 1,  '0.34' ,'5000', 'N1-CN'])
    bondsArray.append([2 ,   3  , 1,  '0.365','5000', 'CN-GLA'])
    bondsArray.append([3 ,   4  , 1,  '0.365','5000', 'GLA-CX'])
    #NB Here the first tail bead need to be saturated
    bondsArray.append([4 ,   5  , 1,  '0.47','5000', 'CX-C1A'])
    bondsArray.append([4 ,   len(tailsArray[0])+index  , 1,  '0.47','5000', 'CX-C1B'])
    anglesArray.append([-1, 'Orient the headgroup and linker'])
    anglesArray.append([1,   2,   3,   1,  '155', '15', 'N1-CN-GLA'])
    anglesArray.append([2,   3,   4,   1,  '130', '20', 'CN-GLA-CX'])
    anglesArray.append([3,   4,   len(tailsArray[0])+index,   1,  '110', '15', 'GLA-CX-C1B'])
    anglesArray.append([3,   4,   5,   1,  '110', '15', 'GLA-CX-C1A'])
    anglesArray.append([5,   4,   len(tailsArray[0])+index,   1,  '125', '25', 'C1A-CX-C1B'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])
    if tailsArray[0][:2] == 'CC':
        anglesArray.append([4,       index,       1+index,   2,     '155',      '25', 'CX-C1A-C2A']) 
    elif tailsArray[0][:2] == 'CD':
        anglesArray.append([4,       index,       1+index,   2,     '140',      '15', 'CX-C1A-D2A'])
    elif tailsArray[0][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[0][:2] == 'DD':
        print ('Not possible yet')
    if tailsArray[1][:2] == 'CC':
        anglesArray.append([4,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '155',      '25', 'CX-C1B-C2B']) 
    elif tailsArray[1][:2] == 'CD':
        anglesArray.append([4,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '140',      '15', 'CX-C1B-D2B'])
    elif tailsArray[1][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[1][:2] == 'DD':
        print ('Not possible yet')
    exclusionsArray.append([-2, '5 {}'.format(len(tailsArray[0])+index)])
    exclusionsArray.append([-2, '[ pairs ]'])
    exclusionsArray.append([-2, '5 {} 1 4.100000e-01    2.350000e+00'.format(len(tailsArray[0])+index)])


if len(headsArray)>0 and headsArray[0]=='MP':
    beadArray.append([1, 'TQ2p', 1, lipidName,  'NP', 1, 1,    ''])
    beadArray.append([2, 'TC2' , 1, lipidName,  'CN', 2, 0,    ''])
    beadArray.append([3, 'SN4a', 1, lipidName,  'GLA', 3, 0,    ''])
    beadArray.append([4, 'SC2' , 1, lipidName,  'CX', 4, 0,  'Branching bead to tail'])
    index += 4
    lCharge += 1.0 # Keep track of overall lipid charge

    bondsArray.append([-1, 'Headgroup and linker bonds'])
    bondsArray.append([1,    2,   1,  '0.34 ', '5000',   'NP-CN'])
    bondsArray.append([2,    3,   1,  '0.365', '5000',   'CN-GLA'])
    bondsArray.append([3,    4,   1,  '0.365', '5000',   'GLA-CX'])
    #NB Here the first tail bead need to be saturated
    bondsArray.append([4,    5,   1,  '0.47', '5000',   'CX-C1A'])
    bondsArray.append([4,    len(tailsArray[0])+index,   1,  '0.47', '5000',   'CX-C1B'])

    anglesArray.append([-1, 'Orient the headgroup and linker'])
    anglesArray.append([1,   2,   3,   1,  '152', '10', 'NP-CN-GLA'])
    anglesArray.append([2,   3,   4,   1,  '137', '15', 'CN-GLA-CX'])
    anglesArray.append([3,   4,   len(tailsArray[0])+index,   1,  '110', '10', 'GLA-CX-C1B'])
    anglesArray.append([3,   4,   5,   1,  '110', '10', 'GLA-CX-C1A'])
    anglesArray.append([5,   4,   len(tailsArray[0])+index,   1,  '125', '25', 'C1A-CX-C1B'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])
    if tailsArray[0][:2] == 'CC':
        anglesArray.append([4,       index,       1+index,   2,     '155',      '25', 'CX-C1A-C2A']) 
    elif tailsArray[0][:2] == 'CD':
        anglesArray.append([4,       index,       1+index,   2,     '140',      '15', 'CX-C1A-D2A'])
    elif tailsArray[0][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[0][:2] == 'DD':
        print ('Not possible yet')
    if tailsArray[1][:2] == 'CC':
        anglesArray.append([4,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '155',      '25', 'CX-C1B-C2B']) 
    elif tailsArray[1][:2] == 'CD':
        anglesArray.append([4,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '140',      '15', 'CX-C1B-D2B'])
    elif tailsArray[1][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[1][:2] == 'DD':
        print ('Not possible yet')
    exclusionsArray.append([-2, '5 {}'.format(len(tailsArray[0])+index)])
    exclusionsArray.append([-2, '[ pairs ]'])
    exclusionsArray.append([-2, '5 {} 1 4.100000e-01    2.350000e+00'.format(len(tailsArray[0])+index)])



if len(headsArray)>0 and headsArray[0]=='B1':
    beadArray.append([1,  'SN3a',  1,   lipidName,  'N1 ',  1,   0, ''])
    beadArray.append([2,  'TN5r',  1,   lipidName,  'GL1',  2,   0, ''])
    beadArray.append([3,  'TN5r',  1,   lipidName,  'GL2',  3,   0, ''])
    beadArray.append([4,  'SC2 ',  1,   lipidName,  'CX ',  4,   0, 'Branching bead to tail'])
    index += 4
    lCharge += 0.0 # Keep track of overall lipid charge

    bondsArray.append([-1, 'Headgroup and linker bonds'])
    bondsArray.append([1,      2,      1,  '0.35', '5000',  'N1-GL1'])
    bondsArray.append([1,      3,      1,  '0.37', '5000',  'N1-GL2'])
    bondsArray.append([2,      3,      1,  '0.35', '5000',  'GL1-GL2'])
    bondsArray.append([2,      4,      1,  '0.37', '5000',  'GL1-CX'])
    bondsArray.append([4,      5,      1,  '0.47', '5000',  'CX-C1A'])
    bondsArray.append([4,    len(tailsArray[0])+index,   1,  '0.47', '5000',  'CX-C1B'])
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])

    anglesArray.append([-1, 'Orient the headgroup and linker'])
    anglesArray.append([1,      2,      4,      2,   '170',   '35',  'N1-GL1-CX'])
    anglesArray.append([1,      3,      4,      2,   '120',   '35',  'N1-GL2-CX'])
    anglesArray.append([3,      2,      4,      2,   '85 ',   '10',  'GL2-GL1-CX'])
    anglesArray.append([2,      4,      5,      2,   '120',   '15',  'GL1-CX-C1A'])
    anglesArray.append([3,      4,      len(tailsArray[0])+index,      2,   '120',   '15',  'GL2-CX-C1B'])
    anglesArray.append([5,      4,      len(tailsArray[0])+index,      2,   '125',   '25',  'C1A-CX-C1B'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])
    if tailsArray[0][:2] == 'CC':
        anglesArray.append([4,       index,       1+index,   2,     '155',      '25', 'CX-C1A-C2A']) 
    elif tailsArray[0][:2] == 'CD':
        anglesArray.append([4,       index,       1+index,   2,     '140',      '15', 'CX-C1A-D2A'])
    elif tailsArray[0][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[0][:2] == 'DD':
        print ('Not possible yet')
    if tailsArray[1][:2] == 'CC':
        anglesArray.append([4,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '155',      '25', 'CX-C1B-C2B']) 
    elif tailsArray[1][:2] == 'CD':
        anglesArray.append([4,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '140',      '15', 'CX-C1B-D2B'])
    elif tailsArray[1][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[1][:2] == 'DD':
        print ('Not possible yet')
    exclusionsArray.append([-2, '2 3 4   ;GL1-GL2-CX'])
    exclusionsArray.append([-2, '3 4     ;GL2-CX'])
    exclusionsArray.append([-2, '1 2 3   ;N1-CN-GL2'])
    exclusionsArray.append([-2, '5 {}    ;C1A-C1B'.format(len(tailsArray[0])+index)])
    exclusionsArray.append([-2, '[ pairs ]'])
    exclusionsArray.append([-2, '5 {} 1   4.100000e-01    2.350000e+00'.format(len(tailsArray[0])+index)])



if len(headsArray)>0 and headsArray[0]=='BP':
    beadArray.append([1,  'SQ2p',  1,   lipidName,  'NP ',  1,   1, ''])
    beadArray.append([2,  'TN5r',  1,   lipidName,  'GL1',  2,   0, ''])
    beadArray.append([3,  'TN5r',  1,   lipidName,  'GL2',  3,   0, ''])
    beadArray.append([4,  'SC2 ',  1,   lipidName,  'CX ',  4,   0, 'Branching bead to tail'])
    index += 4
    lCharge += 1.0 # Keep track of overall lipid charge
    bondsArray.append([-1, 'Headgroup and linker bonds'])
    bondsArray.append([1,      2,      1,  '0.35', '5000',  'N1-GL1'])
    bondsArray.append([1,      3,      1,  '0.37', '5000',  'N1-GL2'])
    bondsArray.append([2,      3,      1,  '0.35', '5000',  'GL1-GL2'])
    bondsArray.append([2,      4,      1,  '0.37', '5000',  'GL1-CX'])
    bondsArray.append([4,      5,      1,  '0.47', '5000',  'CX-C1A'])
    bondsArray.append([4,    len(tailsArray[0])+index,   1,  '0.47', '5000',  'CX-C1B'])
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])
    
    anglesArray.append([-1, 'Orient the headgroup and linker'])
    anglesArray.append([1,      2,      4,      2,   '170',   '35',  'N1-GL1-CX'])
    anglesArray.append([1,      3,      4,      2,   '120',   '35',  'N1-GL2-CX'])
    anglesArray.append([3,      2,      4,      2,   '75 ',   '15',  'GL2-GL1-CX'])
    anglesArray.append([2,      4,      5,      2,   '120',   '15',  'GL1-CX-C1A'])
    anglesArray.append([3,      4,      len(tailsArray[0])+index,      2,   '120',   '15',  'GL2-CX-C1B'])
    anglesArray.append([5,      4,      len(tailsArray[0])+index,   1,  '125', '25', 'C1A-CX-C1B'])

    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])
    if tailsArray[0][:2] == 'CC':
        anglesArray.append([4,       index,       1+index,   2,     '155',      '25', 'CX-C1A-C2A']) 
    elif tailsArray[0][:2] == 'CD':
        anglesArray.append([4,       index,       1+index,   2,     '140',      '15', 'CX-C1A-D2A'])
    elif tailsArray[0][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[0][:2] == 'DD':
        print ('Not possible yet')
    if tailsArray[1][:2] == 'CC':
        anglesArray.append([4,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '155',      '25', 'CX-C1B-C2B']) 
    elif tailsArray[1][:2] == 'CD':
        anglesArray.append([4,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '140',      '15', 'CX-C1B-D2B'])
    elif tailsArray[1][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[1][:2] == 'DD':
        print ('Not possible yet')
    exclusionsArray.append([-2, '2 3 4   ;GL1-GL2-CX'])
    exclusionsArray.append([-2, '3 4     ;GL2-CX'])
    exclusionsArray.append([-2, '1 2 3   ;N1-CN-GL2'])
    exclusionsArray.append([-2, '5 {}    ;C1A-C1B'.format(len(tailsArray[0])+index)])
    exclusionsArray.append([-2, '[ pairs ]'])
    exclusionsArray.append([-2, '5 {} 1   4.100000e-01    2.350000e+00'.format(len(tailsArray[0])+index)])



if len(headsArray)>0 and headsArray[0]=='K2':
    beadArray.append([1,  'TN3a',  1,   lipidName,  'N1 ',  1,   0, ''])
    beadArray.append([2,  'TC2 ',  1,   lipidName,  'CN ',  2,   0, ''])
    beadArray.append([3,  'TN5r',  1,   lipidName,  'GL1',  3,   0, ''])
    beadArray.append([4,  'TN5r',  1,   lipidName,  'GL2',  4,   0, ''])
    beadArray.append([5,  'SC2 ',  1,   lipidName,  'CX ',  5,   0, 'Branching bead to tail'])
    index += 5
    lCharge += 0.0 # Keep track of overall lipid charge

    bondsArray.append([-1, 'Headgroup and linker bonds'])
    bondsArray.append([1,    2,   1,  '0.34', '5000',  'N1-CN'])
    bondsArray.append([2,    3,   1,  '0.34', '5000',  'CN-GL1'])
    bondsArray.append([2,    4,   1,  '0.37', '5000',  'CN-GL2'])
    bondsArray.append([3,    4,   1,  '0.34', '5000',  'GL1-GL2'])
    bondsArray.append([3,    5,   1,  '0.34', '5000',  'GL1-CX'])
    bondsArray.append([4,    5,   1,  '0.47', '5000',  'GL2-CX'])
    #NB Here the first tail bead need to be saturated
    bondsArray.append([5,    6,   1,  '0.47', '5000',  'CX-C1A'])
    bondsArray.append([5,    len(tailsArray[0])+index,   1,  '0.47', '5000',  'CX-C1B'])
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])

    anglesArray.append([-1, 'Orient the headgroup and linker'])
    anglesArray.append([1,   2,   3,   2,  '160', '35',  'N1-CN-GL1'])
    anglesArray.append([1,   2,   4,   2,  '150', '35',  'N1-CN-GL2'])
    anglesArray.append([3,   5,   6,   2,  '130', '15',  'GL1-CX-C1A'])
    anglesArray.append([4,   5,   len(tailsArray[0])+index,   2,  '130', '15',  'GL2-CX-C1B'])
    anglesArray.append([2,   4,   5,   2,  ' 85', '50',  'CN-GL2-CX'])
    anglesArray.append([2,   3,   5,   2,  '130', '50',  'CN-GL1-CX'])
    anglesArray.append([5,      4,   len(tailsArray[0])+index,   1,  '125', '25', 'C1A-CX-C1B'])

    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])
    if tailsArray[0][:2] == 'CC':
        anglesArray.append([5,       index,       1+index,   2,     '155',      '25', 'CX-C1A-C2A']) 
    elif tailsArray[0][:2] == 'CD':
        anglesArray.append([5,       index,       1+index,   2,     '140',      '15', 'CX-C1A-D2A'])
    elif tailsArray[0][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[0][:2] == 'DD':
        print ('Not possible yet')
    if tailsArray[1][:2] == 'CC':
        anglesArray.append([5,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '155',      '25', 'CX-C1B-C2B']) 
    elif tailsArray[1][:2] == 'CD':
        anglesArray.append([5,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '140',      '15', 'CX-C1B-D2B'])
    elif tailsArray[1][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[1][:2] == 'DD':
        print ('Not possible yet')
    exclusionsArray.append([-2, '2 4 3 5 ;CN-GL2-GL1-CX'])
    exclusionsArray.append([-2, '3 4 5   ;GL1-GL2-CX'])
    exclusionsArray.append([-2, '4 5     ;GL2-CX'])
    exclusionsArray.append([-2, '1 2 4   ;N1-CN-GL2'])
    exclusionsArray.append([-2, '5 6 {}   ;CX-C1A-C1B'.format(len(tailsArray[0])+index)])
    exclusionsArray.append([-2, '[ pairs ]'])
    exclusionsArray.append([-2, '6 {} 1   4.100000e-01    2.350000e+00'.format(len(tailsArray[0])+index)])


if len(headsArray)>0 and headsArray[0]=='KP':
    beadArray.append([1,  'TQ2p',  1,   lipidName,  'NP ',  1,   1, ''])
    beadArray.append([2,  'TC2 ',  1,   lipidName,  'CN ',  2,   0, ''])
    beadArray.append([3,  'TN5r',  1,   lipidName,  'GL1',  3,   0, ''])
    beadArray.append([4,  'TN5r',  1,   lipidName,  'GL2',  4,   0, ''])
    beadArray.append([5,  'SC2 ',  1,   lipidName,  'CX ',  5,   0, 'Branching bead to tail'])
    index += 5
    lCharge += 1.0 # Keep track of overall lipid charge

    bondsArray.append([-1, 'Headgroup and linker bonds'])
    bondsArray.append([1,    2,   1,  '0.34', '5000',  'NP-CN'])
    bondsArray.append([2,    3,   1,  '0.34', '5000',  'CN-GL1'])
    bondsArray.append([2,    4,   1,  '0.37', '5000',  'CN-GL2'])
    bondsArray.append([3,    4,   1,  '0.34', '5000',  'GL1-GL2'])
    bondsArray.append([3,    5,   1,  '0.34', '5000',  'GL1-CX'])
    bondsArray.append([4,    5,   1,  '0.47', '5000',  'GL2-CX'])
    #NB Here the first tail bead need to be saturated
    bondsArray.append([5,    6,   1,  '0.47', '5000',  'CX-C1A'])
    bondsArray.append([5,    len(tailsArray[0])+index,   1,  '0.47', '5000',  'CX-C1B'])
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])

    anglesArray.append([-1, 'Orient the headgroup and linker'])
    anglesArray.append([1,   2,   3,   2,  '160', '35',  'N1-CN-GL1'])
    anglesArray.append([1,   2,   4,   2,  '150', '35',  'N1-CN-GL2'])
    anglesArray.append([3,   5,   6,   2,  '130', '25',  'GL1-CX-C1A'])
    anglesArray.append([4,   5,   len(tailsArray[0])+index,   2,  '130', '25',  'GL2-CX-C1B'])
    anglesArray.append([2,   4,   5,   2,  ' 90', '50',  'CN-GL2-CX'])
    anglesArray.append([2,   3,   5,   2,  '130', '50',  'CN-GL1-CX'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])
    anglesArray.append([5,      4,   len(tailsArray[0])+index,   1,  '125', '25', 'C1A-CX-C1B'])

    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])
    if tailsArray[0][:2] == 'CC':
        anglesArray.append([5,       index,       1+index,   2,     '155',      '25', 'CX-C1A-C2A']) 
    elif tailsArray[0][:2] == 'CD':
        anglesArray.append([5,       index,       1+index,   2,     '140',      '15', 'CX-C1A-D2A'])
    elif tailsArray[0][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[0][:2] == 'DD':
        print ('Not possible yet')
    if tailsArray[1][:2] == 'CC':
        anglesArray.append([5,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '155',      '25', 'CX-C1B-C2B']) 
    elif tailsArray[1][:2] == 'CD':
        anglesArray.append([5,       len(tailsArray[0])+index,       len(tailsArray[0])+index+1,   1,     '140',      '15', 'CX-C1B-D2B'])
    elif tailsArray[1][:2] == 'DC':
        print ('Not possible yet')
    elif tailsArray[1][:2] == 'DD':
        print ('Not possible yet')
    exclusionsArray.append([-2, '2 4 3 5 ;CN-GL2-GL1-CX'])
    exclusionsArray.append([-2, '3 4 5   ;GL1-GL2-CX'])
    exclusionsArray.append([-2, '4 5     ;GL2-CX'])
    exclusionsArray.append([-2, '1 2 4   ;N1-CN-GL2'])
    exclusionsArray.append([-2, '5 6 {}  ;CX-C1A-C1B'.format(len(tailsArray[0])+index)])
    exclusionsArray.append([-2, '[ pairs ]'])
    exclusionsArray.append([-2, '6 {} 1   4.100000e-01    2.350000e+00'.format(len(tailsArray[0])+index)])


elif len(headsArray)>0 and headsArray[0]=='NE': # Add DM head neutral
    beadArray.append([1,  'TN3a' ,      1,    lipidName,   'N1' ,        1,      0,     '']) 
    beadArray.append([2,  'SC2'  ,      1,    lipidName,   'CX' ,        2,      0,     '']) 
    beadArray.append([3,  'SN3r',       1,    lipidName,   'GL1',        3,      0,     '']) 
    beadArray.append([4,  'SN3r',       1,    lipidName,   'GL2',        4,      0,     '']) 
    index += 4
    lCharge += 0.0 # Keep track of overall lipid charge

    bondsArray.append([-1, 'Headgroup and linker bonds'])
    bondsArray.append([1,      2,      1,  '0.40', '4000',  'N1-CX'])
    bondsArray.append([2,      3,      1,  '0.43', '3000',  'CX-GL1'])
    bondsArray.append([2,      4,      1,  '0.43', '3000',  'CX-GL2'])
    bondsArray.append([4,      5,      1,  '0.37', '5000',  'GL2-C1A '])
    bondsArray.append([3,   len(tailsArray[0])+index, 1, '0.37', '5000',  'GL1-C1B '])


    anglesArray.append([-1, 'Orient the headgroup and linker'])
    anglesArray.append([1,   2,   4,   2,   '130',   '20',  'N1-CX-GL2'])
    anglesArray.append([1,   2,   3,   2,   '130',   '20',  'N1-CX-GL1'])
    anglesArray.append([4,   5,   6,   2,   '150',   '25',  'GL2-C1A-C2A'])
    anglesArray.append([4,   5,   6,   2,   '150',   '25',  'GL2-C1A-C2A'])
    anglesArray.append([3,   len(tailsArray[0])+index,   len(tailsArray[0])+index+1,   2,   '150',   '25',  'GL1-C1B-C2B'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])



elif len(headsArray)>0 and headsArray[0]=='NP': # Add DM head charged
    beadArray.append([1,  'TQ2p',     1,   lipidName,    'NP',       1,      1, ''])
    beadArray.append([2,  'SC2',      1,   lipidName,    'CX',       2,      0, ''])
    beadArray.append([3,  'SN3r',     1,   lipidName,    'GL1',      3,      0, ''])
    beadArray.append([4,  'SN3r',     1,   lipidName,    'GL2',      4,      0, ''])
    index += 4
    lCharge += 1.0 # Keep track of overall lipid charge

    bondsArray.append([-1, 'Headgroup bonds'])
    bondsArray.append([1,       2,       1,   ' 0.40',    '4000',  'NP-CX'])
    bondsArray.append([2,       3,       1,   ' 0.43',    '3000',  'CX-GL1'])
    bondsArray.append([2,       4,       1,   ' 0.43',    '3000',  'CX-GL2'])
    bondsArray.append([4,       5,       1,    '0.37',    '5000',  'GL2-C1A '])
    bondsArray.append([3,   len(tailsArray[0])+index, 1, '0.37', '5000',  'GL1-C1B '])


    anglesArray.append([-1, 'Orient the headgroup and linker'])
    anglesArray.append([1,   2,   4,   2,   '130',   '20',  'N1-CX-GL2'])
    anglesArray.append([1,   2,   3,   2,   '130',   '20',  'N1-CX-GL1'])
    anglesArray.append([3,   2,   4,   2,   '90',    '10',  'GL1-CX-GL2'])
    anglesArray.append([4,   5,   6,   2,   '150',   '25',  'GL2-C1A-C2A'])
    anglesArray.append([3,   len(tailsArray[0])+index,   len(tailsArray[0])+index+1,   2,   '150',   '25',  'GL1-C1B-C2B'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])



elif len(headsArray)>0 and headsArray[0]=='DT': # Add DL head neutral
    beadArray.append([1,  'SN3a',     1,   lipidName,   'N1 ',      1,      0, ''])
    beadArray.append([2,  'SN4a',     1,   lipidName,   'GL1',      2,      0, ''])
    beadArray.append([3,  'SN4a',     1,   lipidName,   'GL2',      3,      0, ''])
    index += 3
    lCharge += 0.0 # Keep track of overall lipid charge
  
    bondsArray.append([-1, 'Headgroup bonds'])
    bondsArray.append([1,      3,      1,  '0.43',         '3000','N1-GL2'])
    bondsArray.append([1,      2,      1,  '0.43',         '3000','N1-GL1'])
    bondsArray.append([3,      4,      1,  '0.47',         '5000','GL2-C1A, C1A has be a saturated bead']) #For now only SN4a-C1 exists. Missing SN4a-C4
    bondsArray.append([2,      len(tailsArray[0])+index,      1,  '0.47',         '5000','GL1-C1B C1B has to be a saturated bead'])
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])
  
    anglesArray.append([-1, 'Orient the headgroup'])
    anglesArray.append([ 1,      3,      len(tailsArray[0])+index,      2,   '140',   '25',  'N1-GL2-C1B'])
    anglesArray.append([ 1,      2,      4,      2,   '140',   '25',  'N1-GL1-C1A'])
    anglesArray.append([ 2,      1,      3,      2,   '45 ',   '35',  'GL1-N1-GL2'])
    anglesArray.append([ 1,      2,      3,      2,   '80 ',   '55',  'N1-GL1-GL2'])
    if tailsArray[0][:2] == 'CC':
        anglesArray.append([2,   4,    5,  2, '180',      '35', 'GL1-C1A-C2A']) 
    elif tailsArray[0][:2] == 'CD':
        anglesArray.append([2,   4,    5,  2, '180',      '35', 'GL1-C1A-D2A']) 
    elif tailsArray[0][:2] == 'DC':
        anglesArray.append([2,   4,    5,  2, '180',      '35', 'GL1-D1A-C2A']) 
    elif tailsArray[0][:2] == 'DD':
        anglesArray.append([2,   4,    5,  2, '100',      '10', 'GL1-D1A-D2A']) 
    if tailsArray[1][:2] == 'CC':
        anglesArray.append([3,   len(tailsArray[0])+index,    len(tailsArray[0])+index+1,  2, '180',      '35', 'GL2-C1B-C2B']) 
    elif tailsArray[1][:2] == 'CD':
        anglesArray.append([3,   len(tailsArray[0])+index,    len(tailsArray[0])+index+1,  2, '180',      '35', 'GL2-C1B-D2B']) 
    elif tailsArray[1][:2] == 'DC':
        anglesArray.append([3,   len(tailsArray[0])+index,    len(tailsArray[0])+index+1,  2, '180',      '35', 'GL2-D1B-C2B']) 
    elif tailsArray[1][:2] == 'DD':
        anglesArray.append([3,   len(tailsArray[0])+index,    len(tailsArray[0])+index+1,  2, '100',      '10', 'GL2-D1B-D2B']) 
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])


elif len(headsArray)>0 and headsArray[0]=='DP': # Add DL head charged
  beadArray.append([1,  'SQ2p',     1,   lipidName,   'NP ',      1,      1, ''])
  beadArray.append([2,  'SN4a',     1,   lipidName,   'GL1',      2,      0, ''])
  beadArray.append([3,  'SN4a',     1,   lipidName,   'GL2',      3,      0, ''])
  index += 3
  lCharge += 1.0 # Keep track of overall lipid charge

  bondsArray.append([-1, 'Headgroup bonds'])
  bondsArray.append([1,      3,      1,  '0.43',         '3000','NP-GL2'])
  bondsArray.append([1,      2,      1,  '0.43',         '3000','NP-GL1'])
  bondsArray.append([3,      4,      1,  '0.47',         '5000','GL2-C1A, C1A has be a saturated bead']) #For now only SN4a-C1 exists. Missing SN4a-C4
  bondsArray.append([2,      len(tailsArray[0])+index,      1,  '0.47',         '5000','GL1-C1B C1B has to be a saturated bead'])
  bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])

  anglesArray.append([-1, 'Orient the headgroup'])
  anglesArray.append([ 1,      3,      len(tailsArray[0])+index,      2,   '90',   '25',  'NP-GL2-C1B'])
  anglesArray.append([ 1,      2,      4,      2,   '90',   '25',  'NP-GL1-C1A'])
  anglesArray.append([ 2,      1,      3,      2,   '45 ',   '35',  'GL1-NP-GL2'])
  if tailsArray[0][:2] == 'CC':
      anglesArray.append([2,   4,    5,  2,   '180',      '35', 'GL1-C1A-C2A']) 
  elif tailsArray[0][:2] == 'CD':
      anglesArray.append([2,   4,    5,  2,   '180',      '35', 'GL1-C1A-D2A']) 
  elif tailsArray[0][:2] == 'DC':
      anglesArray.append([2,   4,    5,  2,   '180',      '35', 'GL1-D1A-C2A']) 
  elif tailsArray[0][:2] == 'DD':
      anglesArray.append([2,   4,    5,  2,   '100',      '10', 'GL1-D1A-D2A']) 
  if tailsArray[1][:2] == 'CC':
      anglesArray.append([3,   len(tailsArray[0])+index,    len(tailsArray[0])+index+1,  2,   '180',      '35', 'GL2-C1B-C2B']) 
  elif tailsArray[1][:2] == 'CD':
      anglesArray.append([3,   len(tailsArray[0])+index,    len(tailsArray[0])+index+1,  2,   '180',      '35', 'GL2-C1B-D2B']) 
  elif tailsArray[1][:2] == 'DC':
      anglesArray.append([3,   len(tailsArray[0])+index,    len(tailsArray[0])+index+1,  2,   '180',      '35', 'GL2-D1B-C2B']) 
  elif tailsArray[1][:2] == 'DD':
      anglesArray.append([3,   len(tailsArray[0])+index,    len(tailsArray[0])+index+1,  2,   '100',      '10', 'GL2-D1B-D2B']) 
  anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])


# End Add heads
headIndex = index - 1  # To know what was the last headbead
# headIndex = index  # To know what was the last headbead


# Add tail beads
indexToLetter = "A B C D E F G H I J K L M N".split()
for tailIndex in range(0, len(tailsArray)):
    cTail = tailsArray[tailIndex]
    print( tailIndex)
    for cTailBeadIndex in range(0, len(cTail)):
        cTailBead = cTail[cTailBeadIndex]
        # print (cTailBead)
        if (cTailBead not in tailMapp.keys()):
            print("Tail definition \"%s\" not recognized" %
                  (cTailBead), file=sys.stderr)
            sys.exit()
        bType = tailMapp[cTailBead][0]
        atomName = tailMapp[cTailBead][1].replace(
            "??", str(cTailBeadIndex+1) + indexToLetter[tailIndex])
        beadArray.append([index, bType, 1, lipidName, atomName, index, 0, ''])
        # lCharge += 0 # Keep track of overall lipid charge, not needed as current tails are all uncharged

        # Add bond between tail beads
        if cTailBeadIndex > 0:  # linker to first tail alreaddy added
           # if cTailBeadIndex == 1 and \
           #    tailMapp[cTail[cTailBeadIndex-1]][0][0] != 'S' and \
           #    tailMapp[cTail[cTailBeadIndex]][0][0] != 'S':
           #     # If first two tail beads are normal sized (not small) add long tail bond between them
           #     print ('Does not exist')
           #     exit
           #     #bondsArray.append([-3, index-1, index, f"b_{tailMapp[cTail[cTailBeadIndex-1]][3]}_{tailMapp[cTail[cTailBeadIndex]][3]}_mid_5long", ''])
           # else:
            if len(tailsArray) >= 2 and (cTailBeadIndex + 1) == len(cTail):
                bondsArray.append(
                    [-3, index-1, index, f"b_{tailMapp[cTail[cTailBeadIndex-1]][3]}_{tailMapp[cTail[cTailBeadIndex]][3]}", ''])
            else:
                # Add normal mid bond
                bondsArray.append(
                    [-3, index-1, index, f"b_{tailMapp[cTail[cTailBeadIndex-1]][3]}_{tailMapp[cTail[cTailBeadIndex]][3]}", ''])


        # Add angles to support the tails except not for the last tail (which can't have an angle)
        #when cTailBeadIndex is zero, the script provide angle betwen tails
        if ((cTailBeadIndex + 1) < len(cTail) and (cTailBeadIndex != 0) ):
            # print (cTailBeadIndex)
            # print (cTail)
            # print (index)
            anglesArray.append( [-3, index - 1, index, index + 1,  f"a_{tailMapp[cTail[cTailBeadIndex-1]][4]}_{tailMapp[cTail[cTailBeadIndex]][4]}_{tailMapp[cTail[cTailBeadIndex+1]][4]}", ''] )
        # end angle stuff

        index += 1
# End tailsArray loop


# Make .itp headder
itpFile = open(itpFileName, "w")
print(';;;;;; Martini lipid topology for ' +
      lipidCommonName + ', generated using:', file=itpFile)
print('; ' + progString, file=itpFile)
print('; WARNING: Lipids topology was generated following the Martini 3.0 guidelines but not all lipid head group and tail combinations', file=itpFile)
print(';          have been tested; used with care and see <M3lipid main ref> for guidance.\n;', file=itpFile)
print('; Description:', file=itpFile)
print(';   ' + lipidDesc.replace('\\n', '\n;  '), file=itpFile)
current = options["-keyw"].value
if current != None and len(current) > 0:
    print(';@Keywords: '+current, file=itpFile)
print('; Parameterization:', file=itpFile)
print(';   ' + lipidParm.replace('\\n', '\n;  '), file=itpFile)
current = options["-refs"].value
if current != None and len(current) > 0:
    print('; Reference(s): ', file=itpFile)
    print(';   ' + current.replace('\\n', '\n;  '), file=itpFile)
current = options["-crea"].value
if current != None and len(current) > 0:
    print('; Created: ' + current, file=itpFile)
current = options["-auth"].value
if current != None and len(current) > 0:
    print('; Author(s): ' + current, file=itpFile)
current = options["-modi"].value
if current != None and len(current) > 0:
    print('; Modified:', file=itpFile)
    print(';   ' + current.replace('\\n', '\n;  '), file=itpFile)
lArea = options["-area"].value
if lArea != None and len(lArea) > 0:
    print('; Reference area per lipid: ' + lArea + ' nm^2', file=itpFile)
current = options["-warn"].value
if current != None and len(current) > 0:
    print('; Warning(s)/Note(s):', file=itpFile)
    print(';   ' + current.replace('\\n', '\n;  '), file=itpFile)
print(';', file=itpFile)
print('; Molecular topology and mapping of indices:', file=itpFile)

# Add INSANE input string
cTail = lipidTail.replace('c', 'C')  # remove generator specific names
cTail = cTail.replace('F', 'D')  # remove generator specific names
current = '@INSANE alhead='+lipidHead+', altail='+cTail+', alname='+lipidName
if lCharge != None:
    current += ', charge='+str(lCharge)
if lArea != None and len(lArea) > 0:
    current += ', area='+lArea
print(';' + current, file=itpFile)

# Add @RESNTEST, test if using x3 bead resnames how to fine the last letter (e.g. POP is it POPC, POPE, POPS etc)
resntest = ""
cutoflen = 3
if len(lipidName) < cutoflen:  # fix test for short lipid names
    cutoflen = len(lipidName)
if len(headsArray) > 0 and headsArray[0] == 'PI':  # Add PI head\
    resntest = lipidName[0:cutoflen]+'=='+lipidName + \
        ' if: atoms[0]==' + beadArray[0][4] + ' and atoms[4]==GL1'
elif len(headsArray) > 0 and headsArray[0] == 'P1':  # Add PIP_1 head
    resntest = lipidName[0:cutoflen]+'=='+lipidName + \
        ' if: atoms[0]==' + beadArray[0][4] + ' and atoms[4]==P1'
elif len(headsArray) > 0 and headsArray[0] == 'P2':  # Add PIP_2 head
    resntest = lipidName[0:cutoflen]+'=='+lipidName + \
        ' if: atoms[0]==' + beadArray[0][4] + ' and atoms[5]==P2'
elif len(headsArray) > 0 and headsArray[0] == 'P3':  # Add PIP_3 head
    resntest = lipidName[0:cutoflen]+'=='+lipidName + \
        ' if: atoms[0]==' + beadArray[0][4] + ' and atoms[6]==P3'
elif len(headsArray) > 0:  # Else build head (works for simple PC, PE, PA, PG, etc heads)
    resntest = lipidName[0:cutoflen]+'=='+lipidName + \
        ' if: atoms[0]==' + beadArray[0][4]
else:  # If -alhead was empty len(headsArray)==0 then no head to add (DAG, CER, TAG etc)
    resntest = lipidName[0:cutoflen]+'=='+lipidName + \
        ' if: atoms[0]==' + beadArray[0][4]
# Add test if using x3 bead resnames how to fine the last letter (e.g. POP is it POPC, POPE, POPS etc)
print(';@RESNTEST '+resntest, file=itpFile)

# Add all beads
sBeads = ""
beadNameDict = {}
for cBead in beadArray:
    if cBead[0] > 0:
        sBeads += cBead[4].strip()+" "
        beadNameDict[cBead[0]] = cBead[4]
print(';@BEADS '+sBeads, file=itpFile)

# List all bonds
sBonds = ""
for cBond in bondsArray:
    if cBond[0] > 0:
        # Regular bond
        sBonds += beadNameDict[cBond[0]].strip()+"-" + \
            beadNameDict[cBond[1]
                         ].strip()+" "
    elif cBond[0] == -3:
        # Named bond then
        sBonds += beadNameDict[cBond[1]].strip()+"-" + \
            beadNameDict[cBond[2]
                         ].strip()+" "
print(';@BONDS '+sBonds, file=itpFile)

print(';', file=itpFile)
print('', file=itpFile)
print('[moleculetype]', file=itpFile)
print('; molname      nrexcl', file=itpFile)
print('  ' + lipidName + '          1', file=itpFile)

# Write beads
print('\n[atoms]', file=itpFile)
print('; id 	type 	resnr 	residu 	atom 	cgnr 	charge    (mass)', file=itpFile)
for cBead in beadArray:
    if cBead[0] > 0:
        if len(cBead) > 8:  # bead with mass and comment (no option of only mass without comment)
            print('  %2i  %4s  %2i  %4s  %3s  %2i \t%s \t%s \t; %s' % (
                cBead[0], cBead[1], cBead[2], cBead[3], cBead[4], cBead[5], cBead[6], cBead[7], cBead[8]), file=itpFile)
        elif len(cBead[7]) > 0:  # also print comment
            print('  %2i  %4s  %2i  %4s  %3s  %2i \t%s \t; %s' % (
                cBead[0], cBead[1], cBead[2], cBead[3], cBead[4], cBead[5], cBead[6], cBead[7]), file=itpFile)
        else:  # no comment no mass
            print('  %2i  %4s  %2i  %4s  %3s  %2i \t%s' % (
                cBead[0], cBead[1], cBead[2], cBead[3], cBead[4], cBead[5], cBead[6]), file=itpFile)
    elif cBead[0] == -1:  # Regular comment
        print('; ' + cBead[1], file=itpFile)
    elif cBead[0] == -2:  # gromacs system line
        print(cBead[1], file=itpFile)

# Write lipid bonds
print('\n[bonds]', file=itpFile)
# @TODO label will not match every usecase.
print(';  i  j 	name 	(using named bondtypes from martini_v3.0_ffbonded.itp)', file=itpFile)
print(';  i  j 	funct 	force.c.', file=itpFile)
for cBond in bondsArray:
    if cBond[0] > 0:  # Regular bond
        if len(cBond[5]) > 0:  # also print comment
            print('  %2i %2i  %2i \t%s \t%s \t; %s' % (
                cBond[0], cBond[1], cBond[2], cBond[3], cBond[4], cBond[5]), file=itpFile)
        else:  # no comment
            print('  %2i %2i  %2i \t%s \t%s' % (
                cBond[0], cBond[1], cBond[2], cBond[3], cBond[4]), file=itpFile)
    elif cBond[0] == -1:  # Regular comment
        print('; ' + cBond[1], file=itpFile)
    elif cBond[0] == -2:  # gromacs system line
        print(cBond[1], file=itpFile)
    elif cBond[0] == -3:  # named value in ref file
        if len(cBond[4]) > 0:  # also print comment
            print('  %2i %2i \t%s \t; %s' %
                  (cBond[1], cBond[2], cBond[3], cBond[4]), file=itpFile)
        else:  # no comment
            print('  %2i %2i \t%s' %
                  (cBond[1], cBond[2], cBond[3]), file=itpFile)

# Write lipid angles
print('\n[angles]', file=itpFile)
# @TODO label will not match every usecase.
print(';  i  j  k 	name 	(using named angletypes from martini_v3.0_ffbonded.itp)', file=itpFile)
print(';  i  j  k 	funct 	angle 	force.c.', file=itpFile)
for cAngle in anglesArray:
    if cAngle[0] > 0:  # Regular angle
        if len(cAngle[6]) > 0:  # also print comment
            print('  %2i %2i %2i  %2i \t%s \t%s \t; %s' % (
                cAngle[0], cAngle[1], cAngle[2], cAngle[3], cAngle[4], cAngle[5], cAngle[6]), file=itpFile)
        else:  # no comment
            print('  %2i %2i %2i  %2i \t%s \t%s' % (
                cAngle[0], cAngle[1], cAngle[2], cAngle[3], cAngle[4], cAngle[5]), file=itpFile)
    elif cAngle[0] == -1:  # Regular comment
        print('; ' + cAngle[1], file=itpFile)
    elif cAngle[0] == -2:  # gromacs system line
        print(cAngle[1], file=itpFile)
    elif cAngle[0] == -3:  # named value in ref file.
        if len(cAngle[5]) > 0:  # also print comment
            print('  %2i %2i %2i \t%s \t; %s' % (
                cAngle[1], cAngle[2], cAngle[3], cAngle[4], cAngle[5]), file=itpFile)
        else:  # no comment
            print('  %2i %2i %2i \t%s' %
                  (cAngle[1], cAngle[2], cAngle[3], cAngle[4]), file=itpFile)


# Write lipid constraints
if len(constraintsArray) > 0:
    print('\n[constraints]', file=itpFile)
    print(';  i  j  k 	funct 	length', file=itpFile)
    for cConstraint in constraintsArray:
        if cConstraint[0] > 0:
            if len(cConstraint[3]) > 0:  # also print comment
                print('  %2i %2i \t1 \t%s \t; %s' % (
                    cConstraint[0], cConstraint[1], cConstraint[2], cConstraint[3]), file=itpFile)
            else:  # no comment
                print('  %2i %2i \t1 \t%s' % (
                    cConstraint[0], cConstraint[1], cConstraint[2]), file=itpFile)
        elif cConstraint[0] == -1:  # Regular comment
            print('; ' + cConstraint[1], file=itpFile)
        elif cConstraint[0] == -2:  # gromacs system line
            print(cConstraint[1], file=itpFile)

# Write lipid exclusions
if len(exclusionsArray) > 0:
    print('\n[exclusions]', file=itpFile)
    print('; i  j  k  ...', file=itpFile)
    for cExclusion in exclusionsArray:
        if cExclusion[0] > 0:
            print('  %2i %2i \t%s' %
                  (cExclusion[0], cExclusion[1], cExclusion[2]), file=itpFile)
        elif cExclusion[0] == -1:  # Regular comment
            print('; ' + cExclusion[1], file=itpFile)
        elif cExclusion[0] == -2:  # gromacs system line
            print(cExclusion[1], file=itpFile)

print('\n\n', file=itpFile)
itpFile.close()
# End lipid-martini-itp
