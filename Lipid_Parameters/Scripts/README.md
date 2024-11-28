## How to generate itps using script
For generating all the itps automatically use:

> python3 generate-from-table-v02.py -table ILs-CG-Martini3-v1.csv -makeItps yes -ofolder ionizable_lipids_auto_itps

It will make the itps of the lipids specified in the ILs-CG-Martini3-v1.csv file

Include the martini_v3.0_ffbonded.itp file for running. This defines many bonds and angles for the tails.

## Naming scheme

The lipid names will have the form of ABCD, 4 captial letters, where the first two AB will represent the tails. 
The second two, CD, will represent the head and linker combinations. 
For the charged version of the lipid, the D spot will always be name P, for protonated. 

Generally the headgroup amine is named N1 and NP for the neutral and charged lipid, respectively.

An overview of all the lipids can be found in the excel sheet.  


## NB! Branching of tails

ILs with branching tails or more than two tails, can still not yet be constructed with the script. 
Some lipids such as Moderna (SM-102), Pzifer (ALC-0315), Lipid 2, and Lipid 10 have been generated manually from the fragments including branching. 


## insane version

The insane.py has many of the ILs incorporated already.

Working on insane-ILs.py to build all

An exampel
> python2 insane.py -alname test -alhead 'X X P' -allink 'None None None None' -altail 'CCCC CCCC CCC CCCC' -pbc cubix -x 5 -y 5 -z 7 -o mem.gro -l test:100
