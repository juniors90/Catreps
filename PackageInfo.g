# ***************************************************************************
#                                   Catreps
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
#
# This file contains package meta data. For additional information on
# the meaning and correct usage of these fields, please consult the
# manual of the "Example" package as well as the comments in its
# PackageInfo.g file.
#
SetPackageInfo( rec(

PackageName := "Catreps",
Subtitle := "The GAP package for handling representations of groups and categories in positive characteristic.",
Version := "1.0.0",
Date := "20/07/2022", # dd/mm/yyyy format
License := "GPL-2.0-or-later",

Persons := [
  rec(
    FirstNames := "Peter",
    LastName := "Webb",
    WWWHome := "https://www-users.cse.umn.edu/~webb/",
    Email := "webb@math.umn.edu",
    IsAuthor := true,
    IsMaintainer := true,
    PostalAddress := "5000",
    Place := "Minneapolis MN 55455",
    Institution := "School of Mathematics - University of Minnesota",
  ),
  rec(
    FirstNames := "Juan David",
    LastName := "Ferreira",
    WWWHome := "https://github.com/juniors90/",
    Email := "juandavid9a0@gmail.com",
    IsAuthor := false,
    IsMaintainer := true,
    PostalAddress := "5000",
    Place := "Cordoba - Argentina",
    Institution := "Center of Research and Studies in Mathematics - National University of Minnesota",
  ),
],

#SourceRepository := rec( Type := "TODO", URL := "URL" ),
#IssueTrackerURL := "TODO",
PackageWWWHome := "http://www.juniors90.com/",
PackageInfoURL := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
README_URL     := Concatenation( ~.PackageWWWHome, "README.md" ),
ArchiveURL     := Concatenation( ~.PackageWWWHome,
                                 "/", ~.PackageName, "-", ~.Version ),

ArchiveFormats := ".tar.gz",

##  Status information. Currently the following cases are recognized:
##    "accepted"      for successfully refereed packages
##    "submitted"     for packages submitted for the refereeing
##    "deposited"     for packages for which the GAP developers agreed
##                    to distribute them with the core GAP system
##    "dev"           for development versions of packages
##    "other"         for all other packages
##
Status := "dev",

AbstractHTML := Concatenation( [
  "The commands in this package allow you to construct and break apart",
  "representations of both groups and categories, finding their indecomposable",
  "summands and submodule structure. The algorithms of the meataxe for group",
  "representations are included and used where appropriate, but the overall",
  "philosophy is a little different from the meataxe. Methods based on taking", "fixed points are widely used." ] ),

PackageDoc := rec(
  BookName  := "Catreps",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0_mj.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := " The GAP package for handling representations of categories in positive characteristic",
),

Dependencies := rec(
  GAP := ">= 4.11",
  NeededOtherPackages := [ [ "GAPDoc", ">= 1.6.1" ] ],
  SuggestedOtherPackages := [ ],
  ExternalConditions := [ ],
),

AvailabilityTest := ReturnTrue,

TestFile := "tst/testall.g",

Keywords := [ "Groups", "Categories", "Representations", "Positive characteristic" ],

));


