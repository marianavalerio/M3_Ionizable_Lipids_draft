#!/usr/bin/env tclsh
# Turn axis off
display resetview
axes location Off
        
# Check if correct number of arguments are provided
if {[llength $argv] != 2} {
    puts "Usage: $argv0 input.pdb output.stl"
        exit 1
        }

        # Get input and output file paths from command-line arguments
        set pdb_file [lindex $argv 0]
        set stl_file [lindex $argv 1]

        # Load the PDB file
        mol new $pdb_file

        # Set QuickSurf representation with specified parameters
        mol modstyle 0 0 QuickSurf 6.0 0.5 4.0 3.0

        # Render to STL file
        render STL $stl_file true


exit