# ***************************************************************************
#                                   Reps
#       Copyright (C) 2005 Peter Webb 
#       Copyright (C) 2006 Peter Webb, Robert Hank, Bryan Simpkins 
#       Copyright (C) 2007 Peter Webb, Dan Christensen, Brad Froehle
#       Copyright (C) 2020 Peter Webb, Moriah Elkin
#       Copyright (C) 2022 Peter Webb, Ferreira Juan David
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#
#  The overall structure of the reps package was designed and most of it
#  written by Peter Webb <webb@math.umn.edu>, who is also the maintainer. 
#  Contributions were made by
#  Dan Christensen, Roland Loetscher, Robert Hank, Bryan Simpkins,
#  Brad Froehle and others.
# ***************************************************************************
# Read("C:/gap-4.11.1/pkg/Catreps-1.0.0/tst/testall.g");
# This file runs package tests. It is also referenced in the package
# metadata in PackageInfo.g.
#
LoadPackage( "Catreps", "0", false);

TestDirectory(DirectoriesPackageLibrary( "Catreps", "tst" ),
  rec(exitGAP := true));

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
