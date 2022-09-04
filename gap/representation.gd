#! Commands that return a representation:
#! CatRep(category, list of matrices, field) returns a
#! representation of the category over the field, where
#! the generators are represented by the matrices in the list.
DeclareGlobalFunction( "CatRep" );

#! YonedaRep(category, object number, field) returns the
#! representation of the category on the subspace of the
#! category algebra spanned by the morphisms whose domain
#! is the specified object. The representation is projective,
#! by Yoneda's lemma.
DeclareGlobalFunction( "YonedaRep" );

#! YonedaDualRep(category, object number, field) returns the
#! representation of the category on the dual of the subspace
#! of the category algebra spanned by the morphisms whose
#! codomain is the specified object. The representation is
#! injective, because its dual is projective (for the opposite
#! algebra), by Yoneda's lemma.
DeclareGlobalFunction( "YonedaDualRep" );

#! DirectSumRep(rep1, rep2) returns the representation that
#! is the direct sum of the representations rep1 and rep2.
#! DirectSumCatRep(rep1, rep2) returns the representation that is the direct
#! sum of the representations rep1 and rep2 of the category rep1.category.
#!
#! Written by Moriah Elkin July 2018.
#!
#! DirectSumCatRep is called by DirectSumRep(rep1,rep2) when rep1 and rep2 are
#! group representations.
DeclareGlobalFunction( "DirectSumCatRep" );

#! TensorProductRep(rep, rep). . returns a representation
#! which is the tensor product of the two representations.

#! TensorProductCatRep(rep,rep) . . . Kronecker product of two representations
#! Called by TensorProductRep(rep,rep) if rep is a category representation
#!
DeclareGlobalFunction( "TensorProductCatRep" );

#! SubmoduleRep(rep, list of lists of vecs) . . returns the
#! representation which gives the action on the submodule
#! spanned at each object by the corresponding list of vectors.
DeclareGlobalFunction( "SubmoduleRep" );

#! SubmoduleRep(rep, list of lists of vecs) . .  returns the representation which gives
#! the action on the submodule spanned at each object by the corresponding
#! list of vectors. Each list of vectors must be a basis. This is not checked.
#!
#! SubmoduleCatRep is called by SubmoduleRep(rep, list of lists of vecs) when rep is a
#! category representation.
DeclareGlobalFunction( "SubmoduleCatRep" );

#! QuotientRep(rep, list of lists vecs). . returns the
#! representation on the quotient module by the submodule
#! spanned by the vectors. Each list of vectors in the list
#! is a basis for the subspace they span of the representation
#! space for the corresponding object.
DeclareGlobalFunction( "QuotientRep" );


#! QuotientRep(rep, basis structure) . . . returns the representation giving the
#! action on the quotient by an invariant subspace.
#!
#! At the moment this does not work if the basis structure is for the zero subspace.
#!
#! QuotientCatRep is called by QuotientRep(rep, basis structure) when rep is a
#! category representation.
DeclareGlobalFunction( "QuotientCatRep" );

#! ZeroCatRep(cat,field) returns the zero representation of
#! the category.

DeclareGlobalFunction( "ZeroCatRep" );

#! ConstantRep(cat,field) returns the constant (or trivial)
#! representation of the category.
DeclareGlobalFunction( "ConstantRep" );

# -----------------------------------------------------------------+
#! FixedPtBases(n) returns a list of lists of bases. The outer list contains a list for
#! SpinFixedPts of the Yoneda Rep of FI(n) at mathematical object 2, at each object in
#! FI(n); and each of those lists contains the corresponding list of bases at each object
#! (taken at the summand that has fixed points).
#!
#! Written by Moriah Elkin Spring 2019.
DeclareGlobalFunction( "FixedPtBases" );

#! TestYRSummand(basis,eval,summandi) takes in a basis, an evaluation
#! of a Yoneda Rep over GF(2), and the index of the summand of that
#! evaluation that the basis is thought to be equivalent to; it returns
#! whether it is in fact equivalent.
#!
#! Written by Moriah Elkin Spring 2019.
DeclareGlobalFunction( "TestYRSummand" );