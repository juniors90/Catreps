#***************************************************************************
#                                   Reps
#       Copyright (C) 2005 Peter Webb 
#       Copyright (C) 2006 Peter Webb, Robert Hank, Bryan Simpkins 
#       Copyright (C) 2007 Peter Webb, Dan Christensen, Brad Froehle
#       Copyright (C) 2020 Peter Webb, Moriah Elkin
#       Copyright (C) 2022 Ferreira Juan David
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#
#  The overall structure of the reps package was designed and most of it
#  written by Peter Webb <webb@math.umn.edu>, who is also the maintainer. 
#  Contributions were made by
#  Dan Christensen, Roland Loetscher, Robert Hank, Bryan Simpkins,
#  Brad Froehle and others.
#***************************************************************************

#! @Chapter Commands that return a representation

#! @Section Efficient summing over groups

#! @Arguments group, matrices, field(optional)

#! @Returns A representation.

#! @Description This function creates a representation provided the group
#! is not the identity group and the representation is not the zero
#! representation.
#!
#! Modification by Craig Corsi and Peter Webb (University of Minnesota)
#! 2016 to allow the optional field argument.
DeclareGlobalFunction( "Rep" );

#! @Arguments rep

#! @Returns Two fields.

#! @Description Creates fields `rep.permgroup` and `rep.isotopermgroup`
#! which are used by routines which are written to work only with permutation
#! groups, typically involving an algorithm which goes down a stabilizer
#! chain.
DeclareGlobalFunction( "MakeIsoToPermGroup" );

#! @Arguments mat, F

#! @Returns a basis for the nullspace of `mat`.

#! @Description a basis for the nullspace of `mat`, defined 
#! over `F`. This works also when the domain or codomain
#! are the zero vector space.
DeclareGlobalFunction( "SafeNullspaceMat" );

#! @Arguments M

#! @Returns a list a basis vectors for the space spanned
#! by the rows of `M`.

#! @Description If `M` is an empty list, so `M` has no rows, or if
#! the elements of `M` are empty lists, so M has no columns, then the
#! empty list is returned.  
DeclareGlobalFunction( "SafeBaseMat" );

#! @Arguments n, F

#! @Returns an $n \times n$ identity matrix over the field `F`.
#! When `n = 0`, returns an empty list.
DeclareGlobalFunction( "SafeIdentityMat" );

#! @Arguments A, B, n

#! @Returns the product of matrices `A` and `B`,
#! where `A` is $k\times m$ and `B` is $m\times n$.



#! SocleNullspaceMat(matrix, dimension, field)
#! returns a list of vectors forming a basis for the nullspace of the matrix,
#! and deals with the situation where the matrix may have no rows, etc.
#! It is used in SocleSeries(rep); SafeNullspaceMat works in other instances.
#! Written by Peter Webb July 2016.
DeclareGlobalFunction( "SocleNullspaceMat" );



#!
#! RepToMatGroup(rep) returns the matrix group which is the image of the
#! representation.
#!
DeclareGlobalFunction( "RepToMatGroup" );

#! @Arguments group, field

#! @Returns A trivial representation.

#! @Description This works even when the group is the identity group.
DeclareGlobalFunction( "TrivialRep" );

#! @Arguments group, field, n

#! @Returns A trivial representation.

#! @Description trivial action on $n$-dimensional vector space.
#!
#! Code written by Dan Christensen, modified by Peter Webb Aug 2007 to work
#! also with the identity group.
DeclareGlobalFunction( "TrivialRepN" );

#! @Arguments group, field

#! @Returns The zero representation.
DeclareGlobalFunction( "ZeroGroupRep" );

#!
#! FixedPoints. . . . . produce a basis for the fixed points of a module
#!
#! Improvement by Dan Christensen of code originally by Peter Webb.
#!

# DeclareGlobalFunction( "FixedPoints" );

#!
#! FixedQuotient(rep) returns a basis for augmentation ideal . module.
#! Handles trivial groups and zero-dimensional reps.
#!
DeclareGlobalFunction( "FixedQuotient" );

#!
#! SubFixedQuotient(rep, list of vecs)
#! returns a basis for augmentation ideal . submodule spanned by the vecs.
#! The span of the vecs must be a submodule, and this is not checked.
#! Handles trivial groups, zero-dimensional reps and empty u.
#!
DeclareGlobalFunction( "SubFixedQuotient" );



#!
#! TensorProductMatrix(mat,mat) . . . Tensor product of two matrices
#!
#! The GAP command KroneckerProduct seems inexplicably slow and in tests on
#! two 60 x 60 matrices takes about twice as long as the following code.
#!

DeclareGlobalFunction( "TensorProductMatrix" );

#! @Arguments rep1, rep2

#! @Returns A representation which is the tensor product of the two representations.
DeclareGlobalFunction( "TensorProductRep" );

#! @Arguments G1, G2

#! @Returns A

#! @Description TensorProductGroupRep is called by TensorProductRep(rep,rep)
#! when rep is a group representation.
DeclareGlobalFunction( "TensorProductGroupRep" );
DeclareGlobalFunction( "OldTensorProductRep" );

#!
#! TensorProductMorphism(M1,M2) Kronecker product of two morphisms
#!
#! Function introduced August 2007 by Dan Christensen.
#!

DeclareGlobalFunction( "TensorProductMorphism" );
DeclareGlobalFunction( "OldTensorProductMorphism" );

#! @Arguments rep, v

#! @Returns the representation that gives
#! the action on an invariant subspace. The list of
#! vectors must be a basis for an invariant subspace.
#! This is not checked.

#! @Description The code was tidied up by Dan Christensen Aug 2007
#! with the use of List and made to work for the identity group by
#! Peter Webb Aug 2007.
#! In Feb 2020 Peter Webb changed `Size(rep.group)` to `rep.dimension`,
#! apparently correcting an error.
DeclareGlobalFunction( "SubmoduleRep" );

#! SubmoduleGroupRep is called by SubmoduleRep(rep, list of vecs) when rep is a
#! group representation.
DeclareGlobalFunction( "SubmoduleGroupRep" );

#! @Arguments rep, v

#! @Returns The representation giving the
#! action on the quotient by an invariant subspace.

#! @Description The code includes a correction made by Roland Loetscher in March 2007 and was made to work for the identity group by Peter Webb in August 2007.
#!
#! QuotientGroupRep is called by QuotientRep(rep, list of vecs) when rep is a
#! group representation.

DeclareGlobalFunction( "QuotientRep" );
DeclareGlobalFunction( "QuotientGroupRep" );

#!
#! SectionRep(rep, list of vectors, list of vectors) returns the representation
#! on the quotient of the submodule spanned by the first list of vectors, by
#! the submodule spanned by the second list of vectors.
#! The two lists should be independent sets, and should be bases for
#! submodules of rep, the second submodule contained in the first.
#! This is not checked.
#! Written by Peter Webb July 2016.
#!

DeclareGlobalFunction( "SectionRep" );

#!
#!
#!PermutationMatrix(perm,d,field)  We specify d so that the permutation is
#! regarded as permuting [1..d]
#!
#!

DeclareGlobalFunction( "PermutationMatrix" );

#!
#!
#! PermToMatrixGroup( permgrp, field ) . . transforms a permutation group 
#! to a group of permutation matrices
#!
#!

#! PermToMatrixGroup( permgrp, field )
#! 
#!         local   matrix, g, p, d, mgens;
#! 
#!         d := LargestMovedPoint(permgrp);
#!         mgens:=[];
#!         for g in GeneratorsOfGroup(permgrp) do
#!                 matrix:=NullMat(d,d,field);
#!                 for p in [1..d] do
#!                         matrix[p][p^g]:=One(field);
#!                 od;
#!                 Add(mgens,matrix);
#!         od;
#!         return(Group(mgens));
#! end );

#! @Arguments permgrp, field

#! @Returns A `Rep` object as a permutation representation.

#! @Description transforms a permutation group to a permutation
#! representation
DeclareGlobalFunction( "PermGroupToRep" );

#!
#! Pretty print for matrices
#!

DeclareGlobalFunction( "DisplayMatrix" );

#!
#! PrintRep(rep) prints a representation nicely.
#!

DeclareGlobalFunction( "PrintRep" );

#!
#! CanRightRep(group,subgroup,list of elements, element)
#! returns the first element in the list which represents the same
#! coset as the element.
#! It does not verify the validity of its input.
#!

DeclareGlobalFunction( "CanRightRep" );

#!
#! RightCosetReps(group,subgroup) gives a list of representatives  of the 
#! right cosets of a group relative to a given subgroup.
#!

DeclareGlobalFunction( "RightCosetReps" );

#! The function that gives the permutation representation of an element g
#! on a list L of right cosets of a subgroup has a built-in definition as
#! Permutation(g,L,OnRightCosets), but the problem is that L has to be a list 
#! whose elements are themselves lists, each containing the elements of
#! a Right Coset (Cosets are not lists).


#!
#! GLG(n,q) is a generalized version of GeneralLinearGroup
#! which does accept the case of dimension 1. 
#!

DeclareGlobalFunction( "GLG" );

#!
#! RepToGHBI(rep) turns a representation into a GroupHomomorphismByImages
#! with the corresponding data.
#!

DeclareGlobalFunction( "RepToGHBI" );

#! InducedRep(group, subgroup, rep of subgroup) computes the induction
#! of a rep from a subgroup to a group.
#!
#! In more detail, InducedRep(G, H, M) is M tensor_{kH} kG.
#! If m_1, ..., m_n is the standard basis for M and g_1, ..., g_r
#! are the coset representatives of the set {Hg} of right cosets,
#! then the basis chosen for the induced representation is
#!
#!   m_1 tensor g_1, m_2 tensor g_1, ..., m_n tensor g_1,
#!   m_1 tensor g_2, m_2 tensor g_2, ..., m_n tensor g_2,
#!   ...
#!   m_1 tensor g_r, m_2 tensor g_r, ..., m_n tensor g_r,
#!
#! in that order. 
#!

DeclareGlobalFunction( "InducedRep" );

#!
#! InducedMorphism(group, subgroup, matrix) computes the induced
#! homomorphism of a map between representations of the subgroup.
#!
#! Code written by Dan Christensen, University of Western Ontario,August 2007
#!

DeclareGlobalFunction( "InducedMorphism" );

#!
#! InducedInclusion(group, subgroup, rep of subgroup) computes the natural
#! homomorphism from rep to RestrictedRep(G, H, InducedRep(G, H, rep)).
#! With the basis conventions chosen here, this is just a matrix of the
#! form [ I | Z ], where I is an identity matrix and Z is a zero matrix.
#!
#!
#! Code written by Dan Christensen, University of Western Ontario,August 2007
#!

DeclareGlobalFunction( "InducedInclusion" );

#!
#! InducedProjection(group, subgroup, rep of subgroup) computes the natural
#! projection from RestrictedRep(G, H, InducedRep(G, H, rep)) to rep.
#! With the basis conventions chosen here, this is just a matrix of the
#! form [ I ], where I is an identity matrix and Z is a zero matrix.
#!      [ Z ]
#!
#! Code written by Dan Christensen, University of Western Ontario,August 2007
#!

DeclareGlobalFunction( "InducedProjection" );

#!
#! IsHom(rep1, rep2, mat) returns true if the matrix mat represents a
#! homomorphism from rep1 to rep2.  Quite useful for testing.
#!
#! Code written by Dan Christensen, University of Western Ontario,August 2007
DeclareGlobalFunction( "IsHom" );



#! @Arguments G, field

#! @Returns the group ring as a representation.
 
#! @Description The basis is given by the elements of $G$, in the order given
#! by `RightCosetReps(G, Group(Identity(G))`, which is not necessarily
#! the order given by `Elements(G)`.
#!
#! Code written by Dan Christensen, University of Western Ontario, August 2007
DeclareGlobalFunction( "RegularRep" );

#! @Arguments group, field, n

#! @Returns the direct sum of n copies
#! of the group ring as a representation.
#!
#! @Description Code written by Dan Christensen, University of
#! Western Ontario, August 2007
DeclareGlobalFunction( "FreeRep" );

#! @Arguments group, perms, field

#! @Returns The representation on the permutation module determined
#! by the permutation representation in which the group generators
#! act as the permutations in the list.
DeclareGlobalFunction( "PermutationRep" );

#!
#! RelativeTrace(group,subgroup,rep) computes the matrix that corresponds to the mapping
#! tr_Q^P(v)=v(sum gi), where Q is a subgroup
#! of P, v is a vector in a P-representation phi, and the gi are a
#! complete list of representatives of right cosets of Q in P.
#!
#! The real use of this function is to apply the mapping to vectors which
#! are fixed by the subgroup, and then the result is fixed by the whole group.
#!
#! The algorithm is not very clever, and if the index of the subgroup is large
#! it would be better to construct a chain of subgroups between the two groups
#! and compute the relative trace as the product of the relative traces between
#! pairs of groups in the chain. Such an algorithm is used in the command
#! NormRep which returns the relative trace from the identity subgroup.
#!

DeclareGlobalFunction( "RelativeTrace" );



#!
#! OldHomBasis(rep1,rep2) . . returns a basis for the space of 
#! module homomorphisms A -> B. The elements of this basis are matrices.
#!
#! This function is not as fast as the code for HomBasis written by 
#! Dan Christensen, which has a more straightforward setup matrix whose
#! nullspace we find.
#!

DeclareGlobalFunction( "OldHomBasis" );

#!
#! HomBasis(M, N) returns a basis for the space of 
#! kG-module homomorphisms M -> N. The elements of this basis are matrices.
#!
#! The algorithm used finds matrices X so that Xg-gX=0 for all group generators
#! g. This sets up a linear algebra problem more efficiently than solving
#! gXg^-1=X, and was observed independently by John Pliam (1993),
#! Michael Smith (1993) and Dan Christensen (2007).
#! This code written by Dan Christensen, August 2007.
#!
#! GroupHomBasis is called by HomBasis(M, N) when M is a
#! group representation.
#!

DeclareGlobalFunction( "HomBasis" );
DeclareGlobalFunction( "GroupHomBasis" );

#!
#! DimHom(rep1,rep2) . . returns the dimension of the space of 
#! module homomorphisms rep1 -> rep2.
#!

DeclareGlobalFunction( "DimHom" );

#!
#! SumOfImages(rep,rep) . . returns a basis for the sum of images of all 
#! module homomorphisms A -> B
#!
#! GroupSumOfImages is called by SumOfImages(rep,rep) when rep is a
#! group representation.
#!
#! Updated by Peter Webb Sep 11, 2007
#!

DeclareGlobalFunction( "SumOfImages" );
DeclareGlobalFunction( "GroupSumOfImages" );

#!
#! DecomposeSubmodule(repn, basis) . . probably returns a list of two bases 
#! for summands of 
#! the submodule of repn spanned by basis, if the module is decomposable, and
#! if the module is indecomposable it returns a list whose only element is the
#! given basis.
#!
#! GroupDecomposeSubmodule is called by DecomposeSubmodule(repn, basis) when repn
#! is a group representation.
#!
#! This code was written by Bryan Simpkins and Robert Hank, University of
#! Minnesota, April 2006. A related approach was taken by Michael Smith in
#! code written in 1993.
#! It was tidied up by Dan Christensen, University of Western Ontario, Aug 2007.
#! In July 2016 Craig Corsi, University of Minnesota, made a change so that
#! decompositions over fields larger than the field of definition of the 
#! representation are found.
#!

DeclareGlobalFunction( "DecomposeSubmodule" );
DeclareGlobalFunction( "GroupDecomposeSubmodule" );

#!
#! DecomposeOnce(rep) . . probably returns a list of two bases for summands
#! of the representation rep if it is decomposable, and returns a list whose
#! only element is the standard basis if it is indecomposable. 
#!
#! This code was written by Bryan Simpkins and Robert Hank, University of
#! Minnesota, April 2006
#!

DeclareGlobalFunction( "DecomposeOnce" );

#!
#! Decompose(rep) . . returns a list of bases of direct summands of rep which
#! are probably indecomposable.
#!
#! DecomposeGroupRep is called by Decompose(rep) when rep
#! is a group representation.
#!
#! This code was written by Bryan Simpkins and Robert Hank, University of
#! Minnesota, April 2006.
#! It was modified by Dan Christensen and Peter Webb August 2007 so that indecomposable
#! representations are only checked once (before they were checked twice),
#! so as to use space more economically, and so that the result is stored
#! in rep.summands.
#!

DeclareGlobalFunction( "Decompose" );
DeclareGlobalFunction( "DecomposeGroupRep" );

#!
#! MatsOfCosetReps(rep) . . returns a list
#!
#! [rep of a proper subgroup, 
#! list of matrices of right coset representatives of the subgroup,
#! list of right coset representatives of the subgroup].
#!
#! The coset representatives are a Schreier
#! transversal for the stabilizer of a point, and the stabilizer subgroup has Schreier
#! generators for the stabilizer. At the same time matrices which
#! represent the elements which arise are stored.
#!
#! This function is called by NormRep, MatrixOfElement and RestrictedRep.
#! 

DeclareGlobalFunction( "MatsOfCosetReps" );

#!
#! NormRep(rep) . . returns the matrix which represents the sum of the group
#!                  elements.
#!
#! NormRep calls MatsOfCosetReps recursively until the subgroup is the identity group, and
#! multiplies the sums of the matrices of the coset representatives. 
#! 

DeclareGlobalFunction( "NormRep" );
DeclareGlobalFunction( "OldNormRep" );

#!
#! MatrixOfElement(rep,group element) returns the matrix which represents
#! the group element.
#!

DeclareGlobalFunction( "MatrixOfElement" );

#!
#! MatricesOfElements(rep,list of group elements) returns the list of matrices 
#! which represent the group elements.
#!

DeclareGlobalFunction( "MatricesOfElements" );

#!
#! RestrictedRep(group, subgroup, rep of group) computes the restriction
#! of a rep from a group to a subgroup.
#!
#! This code has been rewritten September 2007 by Peter Webb.
#! The algorithm calls the function MatricesOfElements which finds the
#! matrices representing the generators of the subgroup by working down a
#! stabilizer chain for the group and at each stage computing matrices
#! which represent the coset representatives of the stabilizer subgroups
#! and also the generators of the stabilizers. This approach is more 
#! economical that expressing each subgroup generator as a word in the 
#! generators of the big group and then evaluating that word on matrices, because
#! in such words there are repeated subwords which get evaluated again and
#! again.
#!
#! The previous implementation of this function appears as OldRestrictedRep.
#!
#!

DeclareGlobalFunction( "RestrictedRep" );
DeclareGlobalFunction( "OldRestrictedRep" );


#!
#! SymmetricPowerRep(rep, n) computes the action of rep on the n-th symmetric power
#! of the vector space V associated with rep.
#! This routine was written by Brad Froehle at the University of Minnesota, January 2007.
#!

DeclareGlobalFunction( "SymmetricPowerRep" );

#!
#! ProjectiveHomBasis(rep1,rep2) returns a list of matrices which form a basis for
#! the space of module homomorphisms from rep1 to rep2 which factor through
#! a projective module. The algorithm computes the image of the norm map
#! applied to the representation Dual(rep1) tensor rep2.
#!

DeclareGlobalFunction( "ProjectiveHomBasis" );
DeclareGlobalFunction( "OldProjectiveHomBasis" );

#!
#! IsProjectiveMorphism(rep1, rep2, f) returns whether or not
#! the morphism f: rep1 --> rep2 factors through a projective.
#! f is assumed to be a kG-module homomorphism.
#!

DeclareGlobalFunction( "IsProjectiveMorphism" );

#!
#! IsProjectiveRep(rep) returns true if the representation is a projective
#! module, and false otherwise. The algorithm restricts the representation
#! to a Sylow p-subgroup and tests whether |G| times the rank of the norm
#! map equals the dimension of the representation.
#!

DeclareGlobalFunction( "IsProjectiveRep" );

#!
#! ProjectiveFreeSummand(rep) returns a summand of rep which is probably 
#! projective free.  The complementary summand is guaranteed to be 
#! projective, so the returned rep is stably isomorphic to the original rep.
#!

DeclareGlobalFunction( "ProjectiveFreeSummand" );

#!
#! ProjectiveDecomposition(rep) returns a list of two bases, for a submodule
#! which is projective, and for a submodule with probably no non-zero projective summands,
#! whose direct sum is the whole representation. 
#!

DeclareGlobalFunction( "ProjectiveDecomposition" );

#!
#! Interface to the Meataxe.
#!
#! In the following routines several of the meataxe commands available in GAP
#! are converted so that they take representations as arguments and return
#! representations, where the meataxe implementation would have returned
#! meataxe modules. Clearly more of the meataxe commands could be implemented
#! than is done below.
#!

#!
#! MeataxeModuleToRep(rep,meataxemodule) converts a meataxe module to a
#! representation. Because the meataxe module does not store the group
#! being represented, this is obtained from a representation called rep
#! whose only property is that rep.group is the required group.
#!


DeclareGlobalFunction( "MeataxeModuleToRep" );

#!
#! RepToMeataxeModule(rep) converts a representation to a meataxe module.
#!

DeclareGlobalFunction( "RepToMeataxeModule" );

#!
#! ProperSubmodule(rep) calls the meataxe command MTX.ProperSubmoduleBasis
#! and returns a basis for a proper submodule, or [] if there is none.
#!

DeclareGlobalFunction( "ProperSubmodule" );
DeclareGlobalFunction( "IsIrreducibleRep" );
DeclareGlobalFunction( "IsAbsolutelyIrreducibleRep" );
DeclareGlobalFunction( "BasesCompositionSeriesRep" );
DeclareGlobalFunction( "CompositionFactorsRep" );

DeclareGlobalFunction( "SocleRep" );

#!
#!
#! BrauerTraceImage(representation, p-subgroup of rep.group) returns a list
#! of vectors that is a basis for the sum of the images of traces from proper
#! subgroups of the p-group.  Written by Peter Webb June 2016.
#!
#!

DeclareGlobalFunction( "BrauerTraceImage" );

#!
#! FixedPointRep(representation, subgroup of rep.group) returns the 
#! representation of the normalizer of the subgroup on the fixed points
#! of the subgroup.  Written by Peter Webb June 2016.
#!

DeclareGlobalFunction( "FixedPointRep" );

#!
#! BrauerRep(representation, p-subgroup of rep.group) returns the 
#! representation of the normalizer of the subgroup on the fixed points
#! of the subgroup modulo the image of traces from proper subgroups of 
#! the p-group.  Written by Peter Webb June 2016.
#!
#!

DeclareGlobalFunction( "BrauerRep" );



#!
#! SocleSeries(rep) returns a list with two entries [bases, reps] where 
#! bases is a list of bases of the higher socles
#! rep = soc^t(rep), soc^(t-1)(rep), ...
#! in DESCENDING order. The last term in the list is the empty list.
#! reps is a list of the representations on the socle quotients
#! rep/soc^(t-1)(rep), soc^(t-1)(rep)/soc^(t-2)(rep). ...
#! all of which are semisimple. The last term is the last nonzero 
#! representation, and so the list is one shorter than bases.
#! Written by Peter Webb July 2016.
#!

DeclareGlobalFunction( "SocleSeries" );



#! @Arguments rep1, rep2

#! @Returns A representation which is the direct sum of the two representations.
DeclareGlobalFunction( "DirectSumRep" );

#! DirectSumGroupRep(rep1,rep2) returns the representation that is the direct sum of 
#! representations rep1 and rep2 for the same group. Written by Peter Webb February 2020.
#!
#! DirectSumGroupRep is called by DirectSumRep(rep1,rep2) when rep1 and rep2 are
#! group representations.
#!

DeclareGlobalFunction( "DirectSumGroupRep" );

#!
#! ProjectiveFreeCore(rep) returns the representation that is the largest direct summand of 
#! rep with not projective summand. It is only guaranteed to work when the group is a p-group
#! in characteristic p. In other cases it may give a correct answer, and if it does not then an error
#! message is returned. Written by Peter Webb February 2020. There is another function 
#! ProjectiveFreeSummand which invokes Decompose, and because of this will be more limited
#!

DeclareGlobalFunction( "ProjectiveFreeCore" );

#!
#! ProjectiveSummand(rep) returns a basis for a maximal projective direct summand of
#! rep. It is only guaranteed to work when the group is a p-group in characteristic p.
#! In other cases it may give a correct answer, and if it does not then an error
#! message is returned. It does not use Decompose. Written by Peter Webb February 2020.
#!

DeclareGlobalFunction( "ProjectiveSummand" );

#!
#! IsIsomorphicSummand(rep1,rep2) only works when rep1 is indecomposable of dimension prime to the field characteristic. 
#! It returns true if rep1 is isomorphic to a direct summand of rep2, false otherwise. It relies on a result of Benson
#! and Carlson. Written by Peter Webb Feb 2020.
#!

DeclareGlobalFunction( "IsIsomorphicSummand" );

#!
#! KernelIntersection(rep,rep) . . returns a basis for the intersection of the kernels of all 
#! module homomorphisms A -> B
#! Created by Peter Webb May 8, 2019. Corrected April 18, 2020.
#! This uses a version of SafeNullspaceMat with two arguments, the second being the field.
#!

DeclareGlobalFunction( "KernelIntersection" );

#!
#! PrincipalIdempotent(group,prime) returns a vector in the representation space of 
#! RegularRep(group, GF(prime)) that is the block idempotent of the principal block. Thus if b is 
#! this idempotent and e denotes the vector in the regular representation corresponding to the
#! identity elements, the vector b.e is returned. Spinning it gives a basis for the principal block. 
#! This idempotent is constructed using a formula of B. Kulshammer Arch. Math 56 (1991) 313-319.
#! Code written by Peter Webb February 2020.
#!

DeclareGlobalFunction( "PrincipalIdempotent" );

# -------------------------------- NEW FEATURES ---------------------------------------- #
DeclareGlobalFunction( "IrreducibleReps" );
DeclareGlobalFunction( "IsRepr" );
DeclareGlobalFunction( "IrreducibleRepsOfGroup" );
DeclareGlobalFunction( "TensorProductReps" );
DeclareGlobalFunction( "ReprForG" );
DeclareGlobalFunction( "DiagonalRep" );