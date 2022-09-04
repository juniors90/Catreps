#  ---------------------------- categories ------------------
DeclareGlobalFunction( "Domains" );
# ----------------------- homomorphisms ----------------- #
#! SafeNullspaceMat(mat, F) returns a basis for the
#! nullspace Of mat, defined over F. This works also
#! when the domain or codomain are the zero vector space.
DeclareGlobalFunction( "SafeNullspaceMat" );

#! ExtractHom(vec, dimvec1, dimvec2)
#! takes a vector and returns it repackaged as a
#! list of dimvec1[i] by dimvec2[i] matrices.
DeclareGlobalFunction( "ExtractHom" );

# ----------------------- groupsandcateogories -----------------
#! EmptyMat( r, s ) returns an r x s matrix in which each entry is [].
DeclareGlobalFunction( "EmptyMat" );

#! SafeIdentityMat(n, F) returns an nxn identity matrix
#! over the field F. When n = 0, returns an empty list.
DeclareGlobalFunction( "SafeIdentityMat" );

# ----------------------- representation -----------------

#! TensorProductMatrix(mat,mat) . . . Tensor product of two matrices
#!
#! The GAP command KroneckerProduct seems inexplicably slow and in tests on
#! two 60 x 60 matrices takes about twice as long as the following code.
DeclareGlobalFunction( "TensorProductMatrix" );

#! SafeBaseMat(M) returns a list a basis vectors for the space spanned
#! by the rows of M.  If M is an empty list, so M has no rows, or if
#! the elements of M are empty lists, so M has no columns, then the
#! empty list is returned.  
DeclareGlobalFunction( "SafeBaseMat" );

#! SpinFixedPts(rep, obj) returns a list of lists of bases (lists of vectors),
#! where list of bases i is the spin of the fixed points of the group representation
#! that is the i'th summand of the evaluation of rep at object obj; item i in the list
#! is an empty list when there are no generators or no fixed points in that summand of
#! the evaluation.
#!
#! Written by Moriah Elkin Spring 2019.

DeclareGlobalFunction( "SpinFixedPts" );

#! BasisVecDims(n) returns a list of lists of lists of dimensions. The outer list contains a
#! list for SpinBasisVec of the Yoneda Rep of FI(n) at mathematical object 2, at each object
#! in FI(n); and each of those lists contains the dimensions of the corresponding list of
#! bases at each object (where dimension list i contains the dimensions of
#! SpinBasisVec(rep,i)) (an empty list if that summand's basis is empty). I.e.,
#! BasisVecDims[i][j] contains the lengths of the bases generated by spinning the first
#! vector in the basis for the jth summand of the evaluation of the Yoneda rep at object i.
#!
#! Written by Moriah Elkin Spring 2019.
DeclareGlobalFunction( "BasisVecDims" );

#! SpinBasisVec(rep, obj) goes through each summand in the evaluation of rep at obj,
#! takes the first vector in its basis, and returns a basis for the subfunctor generated
#! by this vector (an empty list if the basis is empty). It returns a list of lists of
#! bases, where the i'th list of bases corresponds to the i'th summand of the evaluation.
#! It does not check that the summand does not have fixed points.
#!
#! Written by Moriah Elkin Spring 2019.
DeclareGlobalFunction( "SpinBasisVec" );

#!
#! IsDirectSum(summands,sum) takes in a list of bases (lists of vectors), summands, and a
#! basis (list of vectors), sum, and returns true if the direct sum of summands is sum. 
#!
#! Written by Moriah Elkin Spring 2019.
DeclareGlobalFunction( "IsDirectSum" );

#! TestOneFixedPt(rep) tests whether there is at most one fixed point in the summands
#! of the evaluations of rep at each object. It returns true if there is, false if
#! there is not.
#!
#! Written by Moriah Elkin Spring 2019.
DeclareGlobalFunction( "TestOneFixedPt" );

#! SafeFixedPoints(rep) finds the fixed points of a representation rep; if rep
#! has dimension 0, it returns an empty list.
#!
#! Written by Moriah Elkin Spring 2019.
DeclareGlobalFunction( "SafeFixedPoints" );

#! BrauerTraceImage(representation, p-subgroup of rep.group) returns a list
#! of vectors that is a basis for the sum of the images of traces from proper
#! subgroups of the p-group.  Written by Peter Webb June 2016.
DeclareGlobalFunction( "BrauerTraceImage" );

#! BrauerRep(representation, p-subgroup of rep.group) returns the 
#! representation of the normalizer of the subgroup on the fixed points
#! of the subgroup modulo the image of traces from proper subgroups of 
#! the p-group.  Written by Peter Webb June 2016.
DeclareGlobalFunction( "BrauerRep" );

#! AllSpinDims(n, evalObj,outputRec) spins all vectors in the evaluation at math object
#! evalObj of the Yoneda representation of FI(n) at math object 2 with respect to that
#! Yoneda representation. If outputRec is false, returns a sorted list of the dimensions of
#! the bases; if it is true, returns a record with each dimension list corresponding to the 
#! list of vectors that generated it.
#!
#! Written by Moriah Elkin Spring 2019.
DeclareGlobalFunction( "AllSpinDims" );

#! BasisVecBases(n) returns a list of lists of lists of bases. The outer list contains a 
#! list for SpinBasisVec of the Yoneda Rep of FI(n) at mathematical object 2, at each object
#! in FI(n); and each of those lists contains the corresponding list of bases at each object
#! (where list i contains SpinBasisVec(rep,i)) (an empty list if that summand's basis is
#! empty). I.e., BasisVecBases[i][j] contains the list of bases generated by spinning the
#! first vector in the basis for the jth summand of the evaluation of the Yoneda rep at
#! object i.
#!
#! Written by Moriah Elkin Spring 2019.
DeclareGlobalFunction( "BasisVecBases" );

#! FixedPtDims(n) returns a list of lists of lists of dimensions.
#! The outer list contains a list for SpinFixedPts of the Yoneda
#! Rep of FI(n) at mathematical object 2, at each object
#! in FI(n); and each of those lists contains the dimensions
#! of the corresponding list of bases at each object (an empty
#! list if that summand has no fixed points).
#!
#! Written by Moriah Elkin Spring 2019.
DeclareGlobalFunction( "FixedPtDims" );

#! DimSummands(n,obj) returns the dimensions of the summands in
#! FISummandEvalReps(n,obj,GF(2)): each column of the output is
#! the evaluation at a different mathematical object, each row
#! is a different summand of the YonedaRep of FI(n) at obj, and
#! each element in the matrix is a list with the dimensions of
#! the summands of the evaluation in order.
#!
#! Written by Moriah Elkin Spring 2019.
DeclareGlobalFunction( "DimSummands" );

#! SafeDimHom(rep) returns the dimension of the space of module homomorphisms
#! rep1 -> rep2. If rep has dimension 0, it returns 0.
#!
#! Written by Moriah Elkin Spring 2019.
DeclareGlobalFunction( "SafeDimHom" );

#! ProjSummands(n,obj) returns whether the summands in
#! FISummandEvalReps(n,obj,GF(2)) are projective: each
#! column of the output is the evaluation at a different
#! mathematical object, each row is a different summand
#! of the YonedaRep of FI(n) at obj, and each element in
#! the matrix is a list with whether the summands of the
#! evaluation are projective, in order. Representations
#! with no generators display "fail".
#!
#! Written by Moriah Elkin Spring 2019.
DeclareGlobalFunction( "ProjSummands" );

#! RemoveFromBottom(rep,d) takes representations rep and d.
#! It returns a representation with all images of d removed from the bottom
#! of rep.
#!
#! Written by Moriah Elkin Spring 2019.
DeclareGlobalFunction( "RemoveFromBottom" );

#! RadicalSeries(rep) returns a list with two entries [bases, reps] where 
#! bases is a list of bases of the successive powers of the radical
#! rep = rad^0(rep), rad^1(rep), ...
#! in descending order. The last term in the list is the empty list.
#! reps is a list of the representations on the radical quotients
#! rep/rad^1(rep), rad^1(rep)/rad^2(rep). ...
#! all of which are semisimple. The last term is the last nonzero 
#! representation, and so the list is one shorter than bases.
#! Written by Peter Webb July 2016.
DeclareGlobalFunction( "RadicalSeries" );

#!
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
DeclareGlobalFunction( "SocleSeries" );


#! DisplayButterflyDims(sn,subGens) takes in the symmetric group of a certain dimension
#! and the generators for a subgroup corresponding to a partition. It displays the
#! ButterflyFactors matrix (without decomposing), and returns the ButterflyFactors.
#!
#! Written by Moriah Elkin Spring 2019.
DeclareGlobalFunction( "DisplayButterflyDims" );

#! @Arguments G,H,F

#! @Returns The representation which is the permutation module
#! determined by the right cosets of the subgroup in the group.
#!
#! Modification by Craig Corsi (University of Minnesota)
#! 2016 to make sure the field is correctly recorded.
DeclareGlobalFunction( "PermutationRepOnCosets" );

#! @Description `SafeMatrixMult(A, B, n)` returns the
#! product of matrices `A` and `B`, where `A` is $k\times m$ 
#! and `B` is $m\times n$. $n$ is passed in because if `B` has
#! no rows, then we can not determine $n$.
#!
#! * If $k=0$, returns an empty list, i.e. the matrix with $0$ rows and $n$ columns.
#!
#! * If $n=0$, returns the matrix with $k$ rows and $0$ columns.
#!
#! * If $m=0$, then returns the $k\times n$ zero matrix.
DeclareGlobalFunction( "SafeMatrixMult" );

#!
#! ButterflyFactors(rep, descending filtration, descending filtration) 
#! returns a matrix whose entries are the representations that appear as
#! sections in Zassenhaus' Butterfly Lemma.
#! Each descending filtration is a list of bases of submodules of rep, forming
#! a descending chain. The representations in the output are the factors in
#! a common refinement of the two filtrations, fand their position in the 
#! refinement is indicated by their position in the matrix.
#! Written by Peter Webb July 2016
#!

DeclareGlobalFunction( "ButterflyFactors" );

#! ButterflyDimsRep(rep) takes in a representation, displays the ButterflyFactors matrix
#! (without decomposing), and returns the ButterflyFactors.
#!
#! Written by Moriah Elkin Spring 2019.
DeclareGlobalFunction( "ButterflyDimsRep" );

#! npButterflyDimsRep(rep) takes in a representation and returns the ButterflyFactors
#! (without printing anything).
#!
#! Written by Moriah Elkin Spring 2020.
DeclareGlobalFunction( "npButterflyDimsRep" );

#! ExamineButterflyFactors(b,dRec) takes in a ButterflyFactors matrix and a record of
#! possible factors. It displays the original dimension matrix, and then an ordered list of
#! matrices, where the dimensions in each matrix correspond to the dimensions of the
#! homomorphisms between that element in the ButterflyFactors matrix and the factor in the
#! list, or the number of fixed points in that element in the ButterflyFactors matrix.
#!
#! Written by Moriah Elkin Spring 2019.
DeclareGlobalFunction( "ExamineButterflyFactors" );

#! Pretty print for matrices
DeclareGlobalFunction( "DisplayMatrix" );

#! FISummandEvalReps(n,obj) returns a list of the representations of the summands of the
#! evaluations of the summands of the Yoneda representation over GF(2) of FI(n) at the
#! mathematical object obj. There will be 3 levels of lists: the overall list will contain
#! a list for each submodule of the Yoneda representation, each of which will contain a list
#! of submodule representations for each object at which it has been evaluated.
#!
#! Written by Moriah Elkin Spring 2019.
DeclareGlobalFunction( "FISummandEvalReps" );

# ------------------------- vectors ------------------------------------- #

#! OrthogonalComplement(veclist, dim, field)
#! returns a basis for the orthgonal complement of the space spanned by the
#! list of vectors, in a space of dimension dim (put there in case the list
#! of vectors is empty).
DeclareGlobalFunction( "OrthogonalComplement" );

#!
#! GeneratorDomains(cat) returns a l x l matrix, where l is the number of objects
#! in the category cat, and where the i,j entry is a list of the indices of morphisms
#! from object i to object j.
#!
#! Written by Moriah Elkin August 2018.
#!
DeclareGlobalFunction( "GeneratorDomains" );

# ------------------------- homomorphisms -------------------------------- #
#! RemoveFromTop(rep,d) returns the representation that is the common kernel
#! of all homomorphisms from rep to d.
#!
#! Written by Peter Webb Spring 2020
DeclareGlobalFunction( "RemoveFromTop" );

#! KernelIntersection(rep, rep) . . returns a basis for the intersection of the kernels of all 
#! module homomorphisms A -> B
#! Created by Peter Webb May 8, 2019. Corrected April 18, 2020.
#! This uses a version of SafeNullspaceMat with two arguments, the second being the field.
DeclareGlobalFunction( "KernelIntersection" );

DeclareGlobalFunction( "RadicalRep" );

#! RepToMeataxeModule(rep) converts a representation to a meataxe module.
DeclareGlobalFunction( "RepToMeataxeModule" );

#! @Arguments rep

#! @Returns A a representation which is the dual of rep.

#! @Description Updated Aug 2007 by Peter Webb using the
#! function `InverseGenImages`.
DeclareGlobalFunction( "DualRep" );

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
DeclareGlobalFunction( "RestrictedRep" );

#!
#! MatricesOfElements(rep,list of group elements) returns the list of matrices 
#! which represent the group elements.
DeclareGlobalFunction( "MatricesOfElements" );

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
DeclareGlobalFunction( "RelativeTrace" );

#!
#! IsProjectiveRep(rep) returns true if the representation is a projective
#! module, and false otherwise. The algorithm restricts the representation
#! to a Sylow p-subgroup and tests whether |G| times the rank of the norm
#! map equals the dimension of the representation.
DeclareGlobalFunction( "IsProjectiveRep" );

#!
#! SocleNullspaceMat(matrix, dimension, field)
#! returns a list of vectors forming a basis for the nullspace of the matrix,
#! and deals with the situation where the matrix may have no rows, etc.
#! It is used in SocleSeries(rep); SafeNullspaceMat works in other instances.
#! Written by Peter Webb July 2016.
DeclareGlobalFunction( "SocleNullspaceMat" );

#! InverseGenImages(rep) returns the list of matrices which are the inverses
#! of the matrices in rep.genimages. Two algorithms were tried: in
#! OldInverseGenImages, rather than invert the matrices, they
#! are raised to a power Order(g)-1 for each generator g. In fact, it appears
#! on testing this out that GAP is faster at doing x -> x^-1.
DeclareGlobalFunction( "InverseGenImages" );
DeclareGlobalFunction( "OldInverseGenImages" );

#! @Arguments rep

#! @Returns Two fields.

#! @Description Creates fields `rep.permgroup` and `rep.isotopermgroup`
#! which are used by routines which are written to work only with permutation
#! groups, typically involving an algorithm which goes down a stabilizer
#! chain.
DeclareGlobalFunction( "MakeIsoToPermGroup" );

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
DeclareGlobalFunction( "MatsOfCosetReps" );

#!
#! RightCosetReps(group,subgroup) gives a list of representatives  of the 
#! right cosets of a group relative to a given subgroup.
DeclareGlobalFunction( "RightCosetReps" );

#! RepToGHBI(rep) turns a representation into a GroupHomomorphismByImages
#! with the corresponding data.
DeclareGlobalFunction( "RepToGHBI" );

#! NormRep(rep) . . returns the matrix which represents the sum of the group
#!                  elements.
#!
#! NormRep calls MatsOfCosetReps recursively until the subgroup is the identity group, and
#! multiplies the sums of the matrices of the coset representatives.

DeclareGlobalFunction( "NormRep" );
DeclareGlobalFunction( "OldNormRep" );

#!
#! GLG(n,q) is a generalized version of GeneralLinearGroup
#! which does accept the case of dimension 1. 
DeclareGlobalFunction( "GLG" );

#! TensorProductRep(rep,rep) . . . Kronecker product of two representations
#!
#! TensorProductGroupRep is called by TensorProductRep(rep,rep) when rep is a
#! group representation.

DeclareGlobalFunction( "TensorProductRep" );

