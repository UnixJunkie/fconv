
#============================================================================
# fconv; format conversion, manipulation, and feature computaion
#        of molecular data
#
# Copyright (C) 2007, 2008, 2009, 2010, 2011, 2012 Gerd Neudert
#
# fconv is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, see <http:#www.gnu.org/licenses/>.
#----------------------------------------------------------------------------
#
# author:      Gerd Neudert
#              gneudert(place_at_here)web.de
# supervisor:  Prof. Dr. Gerhard Klebe
#              klebe(place_at_here)staff.uni-marburg.de
# Department of Pharmaceutical Chemistry, Philipps-University of Marburg
#
# We acknowledge the Cambridge Crystallographic Data Centre (CCDC) for
# the permission to use CSD-derived data for bond lengths and bond angles,
# as found in 'ELEMENT_DATA_new_GN.cpp'.
# 
# fconv is a program, intended for parsing and manipulating multiple aspects
# and properties of molecular data. Typical tasks are: Conversion and error
# correction of multiple formats such as PDB, MOL2, SDF, DLG and CIF,
# extracting ligands from PDB as MOL2, automatic or ligand-based cavity
# detection, rmsd calculation and clustering, substructure searches,
# functional alignment and structural superposition, building of crystal
# packings, adding hydrogens, calculation of various properties like the
# number of rotatable bonds, molecular weights or vdW-volumes. The atom type
# classification is based on a consistent assignment of internal atom types,
# which are more differentiated compared to e.g. Sybyl atom types. Apart
# from the predefined mapping of these types onto Sybyl types, the user is
# able to assign own mappings by providing modified template files.
#----------------------------------------------------------------------------
#
# The program is developed since 2007, but I decided to make it open source
# just in 2010 as part of the publication. I appreciate the idea behind free
# software, especially in the field of scientific research. Nevertheless,
# at the moment of publication the documentation inside the code demands
# significant improvement. Furthermore, this program is just a side product
# of my PhD thesis. While having a robust understanding of good programming
# practise nowadays, at the beginning of my PhD I was completely new to C++.
# Thus, there are many parts in the code showing an ugly style and not
# regarding preferable programming rules.
# It's also up to you, to help on improving fconv. I would appreciate, if
# people contact me and share their opinions and changes.
# I am currently using gcc 4.5.x on Ubuntu.
# The code compiles without any errors or warnings with options
# '-pedantic -Wall -O2' on Linux, Windows and MacOS.
# Use '-D_LINUX_OS' if compiling on a linux system.
# 
# Gerd Neudert, Marburg 27.09.2010
#----------------------------------------------------------------------------

 -->  Please see INSTALL.txt for informations on how to compile fconv

 -->  If you have downloaded precompiled binaries, there is only a single
      executable file and you can just start using it.

 -->  fconv should be completely self explaining (Just start it).
       Be sure to READ the HELP carefully (fconv -h) as there is the
       possibility to overwrite your files!

 -->  If you encounter any problems please contact me:

                'gneudert(place_at_here)web.de'

      Please also supply any files you have problems with.
      I would be also pleased for any suggestions on how to improve 'fconv'

 -->  Check for new versions on 'www.agklebe.de'

 -->  Next planned:
 ------------------
   -> reworking the complete commandline interface (to get more consistence
       and less confusion)
   -> reworking the complete programming interface
   -> conversion for different contour maps (e.g. acnt, cnt, grd, map)

Have fun,
27.09.2011, Gerd Neudert
