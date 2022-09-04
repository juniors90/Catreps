#  ***************************************************************************
#                       Catreps
#           Copyright (C) 2008 Peter Webb
#       Copyright (C) 2011 Peter Webb, Fan Zhang
#       Copyright (C) 2020 Moriah Elkin
#       Copyright (C) 2022 Ferreira Juan David
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#
#  The overall structure of the catreps package was designed and most if it
#  written by Peter Webb <webb@math.umn.edu>, who is also the maintainer. 
#  Contributions were made by Dan Christensen,
#  Fan Zhang and Moriah Elkin.
#  ***************************************************************************
#
#! EndomorphismGroup(cat, obj) returns the group of the endomorphisms of obj in
#! cat, in permutation form. It assumes every endomorphism is invertible and that the generators of the 
#! endomorphism group appear among the generators of the category.
#!
#! Written by Moriah Elkin July 2018.
DeclareGlobalFunction( "EndomorphismGroup" );

#! MorphismsRep(rep) returns a l x l matrix, where l is the number of objects
#! in the category of rep, and where the i,j entry is a list of the matrices of
#! morphisms from object i to object j.
#!
#! Written by Moriah Elkin (August 2018), based on code for categories written
#! by Peter Webb.
DeclareGlobalFunction( "MorphismsRep" );

#! Evaluation(rep, obj) returns the representation of the endomorphism
#! group of object obj on the value of the representation rep at obj. (In FI, obj is
#! not the mathematical object in FI, but rather the object number in the category).
#!
#! Written by Moriah Elkin Spring 2019.
DeclareGlobalFunction( "Evaluation" );