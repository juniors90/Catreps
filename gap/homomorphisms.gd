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

#! HomBasis(M, N) returns a basis for the space of 
#! kC-module homomorphisms M -> N. Each element of the output list is a
#! list of matrices, which the ith matrix is a map M[i] -> N[i].
#!
#! The algorithm used finds homomorphisms X so that Xg-gX=0 for all category generators
#! g. The code is an elaboration by Peter Webb (October 2008) of code written 
#! for group representations
#! by Dan Christensen in August 2007.
#!
#! CatHomBasis is called by HomBasis(M, N) when M is a
#! category representation.
DeclareGlobalFunction( "HomBasis" );
DeclareGlobalFunction( "CatHomBasis" );

#! DimHom(rep1,rep2) . . returns the dimension of the space of 
#! module homomorphisms rep1 -> rep2.
DeclareGlobalFunction( "DimHom" );
