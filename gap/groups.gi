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

InstallGlobalFunction( Rep, function(arg)
    # change L by images
    local G, images, rep, arrgh, rho, Ggens;
    arrgh:=arg;
    if not Length(arrgh) in [2..4] or not IsGroup(arrgh[1]) then
        Error("for usage, see ?Rep");
    fi;
    G := arrgh[1];
    images := arrgh[2];
    if images = [] then
        Error("Identity group encountered: creating representations is not properly set up for the identity group.");
        return;    
    fi;
    if   Length(arrgh) = 2 then
        rho := GroupHomomorphismByImages(G, Group( images ) );
    elif Length(arrgh) = 3 then
        Ggens := arrgh[3];
        rho := GroupHomomorphismByImages(G, Group( images ), Ggens, images );
    elif Length(arrgh) = 3 then
        Ggens := arrgh[3];
        rho := GroupHomomorphismByImages(G, Group( images ), Ggens, images );
    fi;
    rep:=rec(
            group:=G,
            generatorsofgroup := GeneratorsOfGroup( G ), # new feature
            rho := rho, # new feature
            rho_images := GeneratorsOfGroup( Images( rho ) ), # new feature
            rho_dom := GeneratorsOfGroup( G ), # new feature
            genimages:=images,
            isRepresentation:=true,
            dimension:=Length(images[1]),
            operations:=GroupRepOps,
        );
    return rep;
    if Length(arrgh) < 4 then
        rep.field:=Field(images[1][1]);
    else
        rep.field:=arrgh[4];
    fi;
end );

InstallGlobalFunction( SafeNullspaceMat, function(mat, F)
    if mat = [] then
        return([]);
    fi;
    if mat[1] = [] then
        return( IdentityMat( Length( mat ), F ) );
    fi;
    return( NullspaceMat(mat) );
end );






InstallGlobalFunction( SafeMatrixMult, function(A, B, n)
    if IsEmpty(A) then     # k = 0
        return [];
    elif n=0 then          # n = 0
        return List(A, row->[]);
    elif IsEmpty(B) then   # m = 0
        return NullMat(Length(A), n);
    else
        return A*B;
    fi;
end );

#!
#! RepToMatGroup(rep) returns the matrix group which is the image of the
#! representation.
#!


InstallGlobalFunction( RepToMatGroup, function(phi)
    return(Group(phi.genimages,IdentityMat(phi.dimension,phi.field)));
end );


InstallGlobalFunction( TrivialRep, function(group, field)
    local mats, onemat;
    onemat:=IdentityMat(1, field);
    mats:=List(GeneratorsOfGroup(group), g->onemat);
    return rec(
        group:=group,
        genimages:=mats,
        field:=field,
        dimension:=1,
        isRepresentation:=true,
        operations:=GroupRepOps
        );
end );


InstallGlobalFunction( TrivialRepN, function(group,field,n)
    local mats,mat;
    mat:=IdentityMat(n,field);
    mats:=List(GeneratorsOfGroup(group), g->mat);
    return rec(
        group:=group,
        genimages:=mats,
        field:=field,
        dimension:=1,
        isRepresentation:=true,
        operations:=GroupRepOps
        );
end );


InstallGlobalFunction( ZeroGroupRep, function(group,field)
    return rec(
     group:=group,
     genimages:=List(GeneratorsOfGroup(group), g->[]),
     field:=field,
     dimension:=0,
     isRepresentation:=true,
     operations:=GroupRepOps
     );
end );





#!
#! FixedQuotient(rep) returns a basis for augmentation ideal . module.
#! Handles trivial groups and zero-dimensional reps.
#!


InstallGlobalFunction( FixedQuotient, function(rep)
    local onemat, v;
    if IsBound(rep.fixedQuotient) then
	    return rep.fixedQuotient;
	fi;
    onemat:=SafeIdentityMat(rep.dimension, rep.field);
    v:=Concatenation(List(rep.genimages, g->g - onemat));
    v:=SafeBaseMat(v);
    rep.fixedQuotient:=v;
    return v;
end );


#!
#! SubFixedQuotient(rep, list of vecs)
#! returns a basis for augmentation ideal . submodule spanned by the vecs.
#! The span of the vecs must be a submodule, and this is not checked.
#! Handles trivial groups, zero-dimensional reps and empty u.
#!


InstallGlobalFunction( SubFixedQuotient, function(rep,u)
    local onemat, v;
    onemat:=SafeIdentityMat(rep.dimension, rep.field);
    v:=Concatenation(List(rep.genimages, g->SafeMatrixMult(u, g - onemat, rep.dimension)));
    return SafeBaseMat(v);
end );




#!
#! TensorProductMatrix(mat,mat) . . . Tensor product of two matrices
#!
#! The GAP command KroneckerProduct seems inexplicably slow and in tests on
#! two 60 x 60 matrices takes about twice as long as the following code.
#!


InstallGlobalFunction( TensorProductMatrix, function( A, B )
    local u, v, matrix;
    matrix := [ ];
    for u in A do
        for v in B do
            Add( matrix, Flat( List( u, x -> x * v ) ) );
        od;
    od;
    return( matrix );
end );

InstallGlobalFunction( TensorProductRep, function( rep1, rep2 )
    return rep1.operations.TensorProductRep( rep1, rep2 );
end );

InstallGlobalFunction( TensorProductCatRep, function(g,h)
    local mgens, i;
    mgens:=[];
    for i in [1..Length(g.genimages)] do
    if IsBound(g.genimages[i]) and IsBound(h.genimages[i]) then
    mgens[i]:=TensorProductMatrix(g.genimages[i],h.genimages[i]);
    fi;
    od;
    return CatRep(g.category,mgens, g.field);
end );



InstallGlobalFunction( OldTensorProductRep, function(g,h)
    local mgens;
    if g.group <> h.group then
        Error("You must have two representations of the same group");
        return;
    fi;
    mgens := List( [1..Length( g.genimages ) ],
                   i -> KroneckerProduct( g.genimages[i], h.genimages[i] ) );
    return Rep(g.group,mgens);
end );


#!
#! TensorProductMorphism(M1,M2) Kronecker product of two morphisms
#!
#! Function introduced August 2007 by Dan Christensen.
#!


InstallGlobalFunction( TensorProductMorphism, function(M1,M2)
    return TensorProductMatrix(M1,M2);
end );

InstallGlobalFunction( OldTensorProductMorphism, function(M1,M2)
    return KroneckerProduct(M1,M2);
end );


#! SubmoduleGroupRep is called by SubmoduleRep(rep, list of vecs) when rep is a
#! group representation.

InstallGlobalFunction( SubmoduleRep, function(rep,v)
    return rep.operations.SubmoduleRep(rep,v);
end );

InstallGlobalFunction( SubmoduleGroupRep, function(rep,v)
    local vs,base,newimages,g;
    vs:=VectorSpace(rep.field,v,[1..rep.dimension]*Zero(rep.field));
    base:=Basis(vs,v);
    newimages:=[];
    for g in rep.genimages do
        Add(newimages, List(base, b->Coefficients(base, b*g)));
    od;
    return rec(
        group:=rep.group,
        genimages:=newimages,
        field:=rep.field,
        dimension:=Length(base),
        isRepresentation:=true,
	operations:=GroupRepOps
	);
end );

#!
#! QuotientGroupRep is called by QuotientRep(rep, list of vecs) when rep is a
#! group representation.
#!
#!

InstallGlobalFunction( QuotientRep, function(rep,v)
    return rep.operations.QuotientRep(rep,v);
end );

InstallGlobalFunction( QuotientGroupRep, function(rep,v)
    local base,d,n,zero,onemat,i,positions,b,p,g,
        mat,newb,newimages,baseinverse, vs, tempbase;
    if Length(v)=0 then
        return rep;
    fi;
    base:=ShallowCopy(BaseMat(v));
    TriangulizeMat(base);
    d:=Length(base);
    n:=rep.dimension;
    if d=n then
        return ZeroGroupRep(rep.group,rep.field);
    fi;
    zero:=Zero(rep.field);
    onemat:=IdentityMat(n,rep.field);
    i:=1;
    positions:=[];
    for b in base do
        while b[i]=zero do
            Add(positions,i);
            i:=i+1;
        od;
        i:=i+1;
    od;
    Append(positions,[i..n]);
    for p in positions do
        Add(base, onemat[p]);
    od;
    baseinverse:=base^-1;
    newimages:=[];
    for g in rep.genimages do
        mat:=[];
        for p in positions do
            b:=g[p]*baseinverse;
            newb:=[];
            for i in [d+1..n] do
                Add(newb,b[i]);
            od;
            Add(mat,newb);
        od;
        Add(newimages, mat);
    od;
    return rec(
        group:=rep.group,
        genimages:=newimages,
        field:=rep.field,
        dimension:=n-d,
        isRepresentation:=true,
        operations:=GroupRepOps
        );
end );

#!
#! SectionRep(rep, list of vectors, list of vectors) returns the representation
#! on the quotient of the submodule spanned by the first list of vectors, by
#! the submodule spanned by the second list of vectors.
#! The two lists should be independent sets, and should be bases for
#! submodules of rep, the second submodule contained in the first.
#! This is not checked.
#! Written by Peter Webb July 2016.
#!
InstallGlobalFunction( SectionRep, function(rep,vectorsa,vectorsb)
    local repa, v, b, newvecsb;
    repa:=SubmoduleRep(rep,vectorsa);
    v:=VectorSpace(rep.field,vectorsa, [1..rep.dimension]*Zero(rep.field));
    b:=Basis(v,vectorsa);
    newvecsb:=List(vectorsb,x->Coefficients(b, x));
    return(QuotientRep(repa,newvecsb));
end );


#!
#!
#!PermutationMatrix(perm,d,field)  We specify d so that the permutation is
#! regarded as permuting [1..d]
#!
#!

InstallGlobalFunction( PermutationMatrix, function(perm,d, field)
    local m, p;
    m:=NullMat(d,d,field);
    for p in [1..d] do
        m[p][p^perm]:=One(field);
    od;
    return(m);
end );


#!
#!
#! PermToMatrixGroup( permgrp, field ) . . transforms a permutation group 
#! to a group of permutation matrices
#!
#!

#! PermToMatrixGroup , function( permgrp, field )
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

#!
#!
#! PermGroupToRep . . transforms a permutation group to a permutation
#!                                representation
#!

InstallGlobalFunction( PermGroupToRep , function( permgrp, field )
        local   matrix, g, p, d, mgens;
        d := LargestMovedPoint(permgrp);
        mgens:=[];
        for g in GeneratorsOfGroup(permgrp) do
            matrix:=NullMat(d,d,field);
            for p in [1..d] do
                matrix[p][p^g] := One(field);
            od;
            Add(mgens,matrix);
        od;
        return(Rep(permgrp, mgens));
end );



#!
#! PrintRep(rep) prints a representation nicely.
#!

InstallGlobalFunction( PrintRep, function( M )
    local g;
    if IsBound( M.name )  then
        Print( M.name );
    elif IsBound(M.longprint) and M.longprint then
        Print( "Representation( ", M.group, ", ", M.genimages, " )\n" );
    else
        Print( "Representation( ", M.group, ", Images \n");
        for g in M.genimages do
            DisplayMatrix(g);
        od;
        Print( " )\n" );
    fi;
end );


#!
#! CanRightRep(group,subgroup,list of elements, element)
#! returns the first element in the list which represents the same
#! coset as the element.
#! It does not verify the validity of its input.
#!


InstallGlobalFunction( CanRightRep, function(G,H,L,g)
   local h,k;
   h:=g^-1;
   for k in L do
     if k*h in H then return k;
     fi; 
   od;
end );

#! The function that gives the permutation representation of an element g
#! on a list L of right cosets of a subgroup has a built-in definition as
#! Permutation(g,L,OnRightCosets), but the problem is that L has to be a list 
#! whose elements are themselves lists, each containing the elements of
#! a Right Coset (Cosets are not lists).






#!
#!
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


InstallGlobalFunction( InducedRep, function(G,H,phi)
    local L,n,F,R,Q,g,i,j,h,S,k,l,ghbiphi;
    ghbiphi:=RepToGHBI(phi);
    R:=[]; n:=phi.dimension; F:=phi.field;
    L:=RightCosetReps(G,H);
    for g in GeneratorsOfGroup(G) do
        Q:=NullMat(n*Length(L),n*Length(L),F);
        for i in [1..Length(L)] do
            j:=Position(L,CanRightRep(G,H,L,L[i]*g));
            h:=L[i]*g*(L[j]^-1);
            S:=ImagesRepresentative(ghbiphi,h);
            for k in [1..n] do
                for l in [1..n] do
                    Q[k+(i-1)*n][l+(j-1)*n]:=S[k][l];
                od;
            od;
        od;
    Add(R, Q);
    od;             
   return Rep(G,R);
end ); 


#!
#! InducedMorphism(group, subgroup, matrix) computes the induced
#! homomorphism of a map between representations of the subgroup.
#!
#! Code written by Dan Christensen, University of Western Ontario,August 2007
#!


InstallGlobalFunction( InducedMorphism, function(G,H,S)
   local index;
   index := Index(G,H);
   return BlockMatrix(List([1..index], i->[i,i,S]), index, index);
   # Alternate implementation:
   # K := figure out the field
   # return KroneckerProduct(IdentityMat(index, K), S)
end );


#!
#! InducedInclusion(group, subgroup, rep of subgroup) computes the natural
#! homomorphism from rep to RestrictedRep(G, H, InducedRep(G, H, rep)).
#! With the basis conventions chosen here, this is just a matrix of the
#! form [ I | Z ], where I is an identity matrix and Z is a zero matrix.
#!
#!
#! Code written by Dan Christensen, University of Western Ontario,August 2007
#!


InstallGlobalFunction( InducedInclusion, function(G,H,phi)
   local dim, index, M, F, i;
   F := phi.field;
   dim := phi.dimension;
   index := Index(G,H);
   M := NullMat(dim, dim*index, F);
   for i in [1..dim] do
      M[i][i] := One(F);
   od;
   return M;
end );


#!
#! InducedProjection(group, subgroup, rep of subgroup) computes the natural
#! projection from RestrictedRep(G, H, InducedRep(G, H, rep)) to rep.
#! With the basis conventions chosen here, this is just a matrix of the
#! form [ I ], where I is an identity matrix and Z is a zero matrix.
#!      [ Z ]
#!
#! Code written by Dan Christensen, University of Western Ontario,August 2007
#!


InstallGlobalFunction( InducedProjection, function(G,H,phi)
   local dim, index, M, F, i;
   F := phi.field;
   dim := phi.dimension;
   index := Index(G,H);
   M := NullMat(dim*index, dim, F);
   for i in [1..dim] do
      M[i][i] := One(F);
   od;
   return M;
end );


#!
#! IsHom(rep1, rep2, mat) returns true if the matrix mat represents a
#! homomorphism from rep1 to rep2.  Quite useful for testing.
#!
#! Code written by Dan Christensen, University of Western Ontario,August 2007
#!


InstallGlobalFunction( IsHom , function(rep1, rep2, mat)
    local i;
    if rep1.group <> rep2.group then
        Error("You must have two representations of the same group");
        return;
    fi;
    for i in [1..Length(rep1.genimages)] do
        if rep1.genimages[i] * mat <> mat * rep2.genimages[i] then
            return false;
        fi;
    od;
    return true;
end );

InstallGlobalFunction( PermutationRepOnCosets, function(G,H,F)
   local n,L,M,phi,g,i,j,R,glg,k;
   L:=[]; R:=[];
   L:=RightCosets(G,H);
   n:=Length(L);
   for k in [1..Length(GeneratorsOfGroup(G))] do
     g:=GeneratorsOfGroup(G)[k];
     M:=NullMat(n,n,F);
     for i in [1..n] do
       j:=0;
       repeat j:=j+1;
       until L[i]*g=L[j];
       M[i][j]:=Z(Size(F))^0;
     od;
     R[k]:=M;
   od;
   phi:=Rep(G,R,F);
   return phi;
end );


InstallGlobalFunction( RegularRep, function(G, field)
    return PermutationRepOnCosets(G, Group(Identity(G)), field);
end );


InstallGlobalFunction( FreeRep, function(G, field, n)
    return TensorProductRep(TrivialRepN(G, field, n), RegularRep(G, field));
end );

InstallGlobalFunction( PermutationRep, function(group,perms,field)
    local d;
    d:=Maximum(List(perms,LargestMovedPointPerm));
    return Rep(group,List(perms, x->PermutationMatrix(x,d,field)));
end );

#!
#! OldHomBasis(rep1,rep2) . . returns a basis for the space of 
#! module homomorphisms A -> B. The elements of this basis are matrices.
#!
#! This function is not as fast as the code for HomBasis written by 
#! Dan Christensen, which has a more straightforward setup matrix whose
#! nullspace we find.
#!


InstallGlobalFunction( OldHomBasis, function(g,h)
local  f, basis, i, j, m, u;
f:=FixedPoints@(TensorProductRep(DualRep(g),h));
basis:=[];
u:=[];
for m in f do
    for i in [1..g.dimension] do
    u[i]:=[];
        for j in [1..h.dimension] do
            u[i][j]:=m[h.dimension*(i-1)+j];
        od;
    od;
    Add(basis,ShallowCopy(u));
od;
return(basis);
end );


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


InstallGlobalFunction( HomBasis, function(M,N)
    return M.operations.HomBasis(M,N);
end );

InstallGlobalFunction( GroupHomBasis, function(M,N)
    local dimM, dimN, r, i, j, k, l, v, f, basis, m, u;
    dimM := M.dimension;
    dimN := N.dimension;
    r := Length(M.genimages);
    v:=NullMat(dimM*dimN, dimM*dimN*r, M.field);
    for i in [1..dimM] do
        for k in [1..dimN] do
            for l in [1..r] do
                for j in [1..dimM] do
                    v[(j-1)*dimN+k][(l-1)*dimM*dimN+(i-1)*dimN+k] :=
                    v[(j-1)*dimN+k][(l-1)*dimM*dimN+(i-1)*dimN+k] + M.genimages[l][i][j];
                od;
                for j in [1..dimN] do
                    v[(i-1)*dimN+j][(l-1)*dimM*dimN+(i-1)*dimN+k] :=
                    v[(i-1)*dimN+j][(l-1)*dimM*dimN+(i-1)*dimN+k] - N.genimages[l][j][k];
                od;
            od;
        od;
    od;
    f := NullspaceMat(v);
    basis:=[];
    for m in f do
        u:=[];
        for i in [1..M.dimension] do
            u[i]:=[];
            for j in [1..N.dimension] do
                u[i][j]:=m[N.dimension*(i-1)+j];
            od;
        od;
        Add(basis, u);
    od;
    return(basis);
end );


#!
#! DimHom(rep1,rep2) . . returns the dimension of the space of 
#! module homomorphisms rep1 -> rep2.
#!


InstallGlobalFunction( DimHom, function(g,h)
    return Length(HomBasis(g,h));
end );


InstallGlobalFunction( GroupSumOfImages, function(g,h)
    return SafeBaseMat(Concatenation(HomBasis(g,h)));
end );





InstallGlobalFunction( GroupDecomposeSubmodule, function(repn, basis)
	local newrep, initlist, a, kernel, image, b, d, n, z;
	if Length(basis) <= 1 then return [basis]; fi;
	newrep := SubmoduleRep(repn, basis);
	z:=PrimitiveElement(repn.field);
	initlist := HomBasis(newrep,newrep);
	Add(initlist, 0 * initlist[1]);
	b := Length(initlist);
	while b > 0 do;
		d := b - 1;
		while d > 0 do;
			a := z*initlist[b] + initlist[d];
			n := 1;
			while n < Length(basis) do;
				a:= a*a;
				n:= 2*n;
			od;
			kernel := NullspaceMat(a);
			if not(Length(kernel) = 0 or Length(kernel) = Length(basis)) then
				image  := BaseMat(a);
				return [kernel * basis, image * basis];
			fi;
			d := d - 1;
		od;
		b := b - 1;
	od;
	return [basis];
end );


#!
#! DecomposeOnce(rep) . . probably returns a list of two bases for summands
#! of the representation rep if it is decomposable, and returns a list whose
#! only element is the standard basis if it is indecomposable. 
#!
#! This code was written by Bryan Simpkins and Robert Hank, University of
#! Minnesota, April 2006
#!


InstallGlobalFunction( DecomposeOnce, function(rep)
	return DecomposeSubmodule(rep, IdentityMat(rep.dimension, rep.field));
end );


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


InstallGlobalFunction( Decompose, function(rep)
    return rep.operations.Decompose(rep);
end );

InstallGlobalFunction( DecomposeGroupRep, function(rep)
	local summands, result, q;
	if IsBound(rep.summands) then return(rep.summands);
	fi;
	summands := [IdentityMat(rep.dimension,rep.field)];
	q := 1;
	# We maintain the following invariants:
	# - summands is a list of lists of vectors; the union of these
	#   lists forms a basis for rep.field^(rep.dimension).
        # - the summands at positions < q appear to be indecomposable;
	#   those at positions >= q haven't been investigated completely.
	while IsBound(summands[q]) do;
		result := DecomposeSubmodule(rep, summands[q]);
		if Length(result) = 2 then
			summands[q] := result[1];
			Add(summands, result[2]);
		else
			q := q + 1;
		fi;
	od;
	rep.summands:=summands;
	return summands;
end );


#!
#! MatrixOfElement(rep,group element) returns the matrix which represents
#! the group element.
#!


InstallGlobalFunction( MatrixOfElement, function(rep,g)
local newrep, newg, mat, s, temp;
    MakeIsoToPermGroup(rep);
    newrep:=rec(
        group:=rep.permgroup,
        genimages:=rep.genimages,
        field:=rep.field,
        dimension:=rep.dimension,
        isRepresentation:=true,
        operations:=GroupRepOps
    );
    newg:=Image(rep.isotopermgroup,g);
    mat:=IdentityMat(rep.dimension,rep.field);
    while Size(newrep.group)>1 do
        s:=SmallestMovedPoint(newrep.group);
        temp:=MatsOfCosetReps(newrep);
        newrep:=temp[1];
        mat:=temp[2][s^newg]*mat;
        newg:=newg*temp[3][s^newg]^-1;
    od;
    return(mat);
end );


InstallGlobalFunction( OldRestrictedRep, function(G,H,phi)
   local ghbiphi,R;
   ghbiphi:=RepToGHBI(phi);
   R:=List(GeneratorsOfGroup(H), g->ImagesRepresentative(ghbiphi, g));
    return rec(
        group:=H,
        genimages:=R,
        field:=phi.field,
        dimension:=phi.dimension,
        isRepresentation:=true,
        operations:=GroupRepOps
        );
end );



#!
#! SymmetricPowerRep(rep, n) computes the action of rep on the n-th symmetric power
#! of the vector space V associated with rep.
#! This routine was written by Brad Froehle at the University of Minnesota, January 2007.
#!


InstallGlobalFunction( SymmetricPowerRep , function(rep, n)
	local Sn;

	#! Define a function Sn which can be called repeatedly for each matrix corresponding to a generator of the group.
	Sn := function(A, n)
		local MyProduct, DualList, dimV, expEnd, dimOut, output, waysToSum, dual, i, j, l;

		#! MyProduct(A, x, y)
		#! 
		#! Assumptions: A is a square n-by-n matrix, x and y are lists of the same length
		#! whose entries are elements of [1..n]
		#!
		#! Output: the product A[x[1]][y[1]] * A[x[2]][y[2]] * A[x[3]][y[3]] * ...
		#! 
		MyProduct := function(A, x, y)
			local output, i;
			output := 1;
			for i in [1..Length(x)]
				do
				output := output * A[x[i]][y[i]];
			od;
			return output;
		end;
		
		#! DualList(A) 
		#!
		#! Assumptions: A is a list.
		#!
		#! Output: A list in which 1 is repeated A[1] times, 2 is repeated A[2] times, etc...
		#!
		DualList := function(A)
			local i, output;
			output := [];
			for i in [1..Length(A)]
				do
				Append(output, List([1..A[i]], x->i));
			od;
			return output;
		end;

		#! Define some variables which are used repeatedly:
		#!  * dimV is the dimension of the initial representation.
		#!  * expEnd is a list whose elements describe all possible monomials of total degree dimV.
		#!  * dimOut is the dimension of the output representation, i.e. the number of monomials in expEnd.
		#!
		dimV := Size(A);
		expEnd := Reversed(OrderedPartitions(dimV+n,dimV)-1);
		dimOut := Size(expEnd);
		
		#! Initialize output to a dimOut x dimOut matrix of all zeros
		#!
		output := [];
		for i in [1..dimOut]
			do
			output[i] := List([1..dimOut],x->0);
		od;

		#! Iterate, calculating each entry of output individually.
		#!
		for i in [1..dimOut]
			do
			# Because this next call (PermutationsList) might be slow, we calculate it as few times as possible.
			waysToSum := PermutationsList(DualList(expEnd[i]));
			for j in [1..dimOut]
				do
				# Calculate the ji-th entry here
				dual := DualList(expEnd[j]);
				for l in waysToSum
					do
					output[j][i] := output[j][i] + MyProduct(A, dual, l);
				od;
			od;
		od;
		return output;
	end;
	#! End function Sn.
	
	return Rep(rep.group, List(rep.genimages, x-> Sn(x,n)));
end );


#!
#! ProjectiveHomBasis(rep1,rep2) returns a list of matrices which form a basis for
#! the space of module homomorphisms from rep1 to rep2 which factor through
#! a projective module. The algorithm computes the image of the norm map
#! applied to the representation Dual(rep1) tensor rep2.
#!


InstallGlobalFunction( ProjectiveHomBasis, function( rep1, rep2 )
    local f, basis, i, j, m, u;
    f:=SafeBaseMat( NormRep( TensorProductRep( DualRep( rep1 ), rep2 ) ) );
    basis:=[];
    for m in f do
        u:=[];
        for i in [1..rep1.dimension] do
        u[i]:=[];
            for j in [1..rep2.dimension] do
                u[i][j]:=m[rep2.dimension*(i-1)+j];
            od;
        od;
        Add(basis,u);
    od;
    return(basis);
end );

InstallGlobalFunction( OldProjectiveHomBasis, function(rep1,rep2)
    local  f, basis, i, j, m, u;
    f:=SafeBaseMat(OldNormRep(TensorProductRep(DualRep(rep1),rep2)));
    basis:=[];
    u:=[];
    for m in f do
        for i in [1..rep1.dimension] do
        u[i]:=[];
            for j in [1..rep2.dimension] do
                u[i][j]:=m[rep2.dimension*(i-1)+j];
            od;
        od;
        Add(basis,ShallowCopy(u));
    od;
    return(basis);
end );



#!
#! IsProjectiveMorphism(rep1, rep2, f) returns whether or not
#! the morphism f: rep1 --> rep2 factors through a projective.
#! f is assumed to be a kG-module homomorphism.
#!


InstallGlobalFunction( IsProjectiveMorphism, function(M,N,f)
    local mat, n;
    mat:=ProjectiveHomBasis(M,N);
    n:=Length(mat);
    Add(mat,f);
    mat:=List(mat,Flat);
    return RankMat(mat)=n;
end );


#!
#! ProjectiveFreeSummand(rep) returns a summand of rep which is probably 
#! projective free.  The complementary summand is guaranteed to be 
#! projective, so the returned rep is stably isomorphic to the original rep.
#!


InstallGlobalFunction( ProjectiveFreeSummand, function(M)
    local comps, comp, pfbasis;
    comps:=Decompose(M);
    pfbasis:=[];
    for comp in comps do
        if not IsProjectiveRep(SubmoduleRep(M, comp)) then
            Append(pfbasis, comp);
        fi;
    od;
    return SubmoduleRep(M, pfbasis);
end );


#!
#! ProjectiveDecomposition(rep) returns a list of two bases, for a submodule
#! which is projective, and for a submodule with probably no non-zero projective summands,
#! whose direct sum is the whole representation. 
#!


InstallGlobalFunction( ProjectiveDecomposition, function(M)
    local comps, comp, pfbasis, pbasis;
    comps:=Decompose(M);
    pfbasis:=[];
    pbasis:=[];
    for comp in comps do
        if not IsProjectiveRep(SubmoduleRep(M, comp)) then
            Append(pfbasis, comp);
        else Append(pbasis, comp);
        fi;
    od;
    return ([pbasis,pfbasis]);
end );

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


InstallGlobalFunction( MeataxeModuleToRep, function(rep,meataxemodule)
    return(rec(group:=rep.group,
           genimages:=MTX.Generators(meataxemodule),
           field:=MTX.Field(meataxemodule),
           dimension:=MTX.Dimension(meataxemodule),
           isRepresentation:=true,
           operations:=GroupRepOps)
        );
end );





#!
#! ProperSubmodule(rep) calls the meataxe command MTX.ProperSubmoduleBasis
#! and returns a basis for a proper submodule, or [] if there is none.
#!


InstallGlobalFunction( ProperSubmodule, function(rep)
    local basis;
    if rep.dimension=0 then return([]); fi;
    basis:=(MTX.ProperSubmoduleBasis(RepToMeataxeModule(rep)));
    if basis=fail then return([]); fi;
    return(basis);
end );

InstallGlobalFunction( IsIrreducibleRep, function(rep)
    if rep.dimension=0 then return(false); fi;
    return(MTX.IsIrreducible(RepToMeataxeModule(rep)));
end );

InstallGlobalFunction( IsAbsolutelyIrreducibleRep, function(rep)
    if rep.dimension=0 then return(false); fi;
    return(MTX.IsAbsolutelyIrreducible(RepToMeataxeModule(rep)));
end );

InstallGlobalFunction( BasesCompositionSeriesRep, function(rep)
    if rep.dimension=0 then return([[]]); fi;
    return(MTX.BasesCompositionSeries((RepToMeataxeModule(rep))));
end );

InstallGlobalFunction( CompositionFactorsRep, function(rep)
    local modules;
    if rep.dimension=0 then return([]); fi;
    modules:=MTX.CompositionFactors(RepToMeataxeModule(rep));
    return(List(modules,x->MeataxeModuleToRep(rep,x)));    
end );



InstallGlobalFunction( SocleRep, function(rep)
    if rep.dimension=0 then return([]); fi;
    return(MTX.BasisSocle(RepToMeataxeModule(rep)));
end );



#!
#!
#! FixedPointRep(representation, subgroup of rep.group) returns the 
#! representation of the normalizer of the subgroup on the fixed points
#! of the subgroup.  Written by Peter Webb June 2016.
#!
#!

InstallGlobalFunction( FixedPointRep, function(rep,gp)
    local n, resn, resgp;
    n:=Normalizer(rep.group,gp);
    resn:=RestrictedRep(rep.group,n,rep);
    resgp:=RestrictedRep(rep.group,gp,rep);
    return(SubmoduleRep(resn,FixedPoints@(resgp)));
end );

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
InstallGlobalFunction( ButterflyFactors, function(rep,L1,L2)
    local factors, i, j, vecsa, vecsb;
    factors:=[];
    for i in [1..Length(L1)-1] do
        factors[i]:=[];
        for j in [1..Length(L2)-1] do
            vecsa:=SumIntersectionMat(L1[i],L2[j])[2];
            vecsb:=SumIntersectionMat(SumIntersectionMat(L1[i+1],L2[j])[2],
            SumIntersectionMat(L1[i],L2[j+1])[2])[1];
            factors[i][j]:=SectionRep(rep,vecsa,vecsb);
        od;
    od;
    return(factors);
end );

InstallGlobalFunction( DirectSumRep, function(rep1,rep2)
    return rep1.operations.DirectSumRep(rep1,rep2);
end );


#!
#! DirectSumGroupRep(rep1,rep2) returns the representation that is the direct sum of 
#! representations rep1 and rep2 for the same group. Written by Peter Webb February 2020.
#!
#! DirectSumGroupRep is called by DirectSumRep(rep1,rep2) when rep1 and rep2 are
#! group representations.
#!
#!


InstallGlobalFunction( DirectSumGroupRep, function( rep1, rep2 )
    local i, gimages, x, zerovec1, zerovec2, padded1, padded2;
    if rep1.field<>rep2.field or rep1.group<>rep2.group then
        Error("Representations incompatible.");
    fi;
    zerovec1:=Zero(rep1.field)*[1..rep1.dimension];
    zerovec2:=Zero(rep2.field)*[1..rep2.dimension];
    gimages:=[];
    for i in [1..Length(rep1.genimages)] do
    padded1:=List(rep1.genimages[i],x->Concatenation(x,zerovec2));
    padded2:=List(rep2.genimages[i],x->Concatenation(zerovec1,x));
    gimages[i]:=Concatenation(padded1,padded2);
    od;
    return rec(
               group:=rep1.group,
               genimages:=gimages,
               field:=rep1.field,
               dimension:=rep1.dimension+rep2.dimension,
               isRepresentation:=true,
               operations:=GroupRepOps
            );
end );


#!
#! ProjectiveFreeCore(rep) returns the representation that is the largest direct summand of 
#! rep with not projective summand. It is only guaranteed to work when the group is a p-group
#! in characteristic p. In other cases it may give a correct answer, and if it does not then an error
#! message is returned. Written by Peter Webb February 2020. There is another function 
#! ProjectiveFreeSummand which invokes Decompose, and because of this will be more limited
#!


InstallGlobalFunction( ProjectiveFreeCore, function(rep)
    local n, resrep, vectors, sub;
    resrep:=RestrictedRep(rep.group, SylowSubgroup(rep.group,Characteristic(rep.field)),rep);
    n:=NormRep(resrep);
    if Length(n)=0 then
        return(rep);
    fi;
    vectors:=BaseSteinitzVectors(IdentityMat(rep.dimension,rep.field),NullspaceMat(n));
    sub:=Spin(rep,vectors.factorspace);
    if IsProjectiveRep(SubmoduleRep(rep,sub))=false then
        Error("The routine tried to factor out a non-projective module. It only works for p-groups.");
        return;
    fi;
    return(QuotientRep(rep,sub));
end );

#!
#! ProjectiveSummand(rep) returns a basis for a maximal projective direct summand of
#! rep. It is only guaranteed to work when the group is a p-group in characteristic p.
#! In other cases it may give a correct answer, and if it does not then an error
#! message is returned. It does not use Decompose. Written by Peter Webb February 2020.
#!

InstallGlobalFunction( ProjectiveSummand, function(rep)
    local n, resrep, vectors, sub;
    resrep:=RestrictedRep(rep.group, SylowSubgroup(rep.group,Characteristic(rep.field)),rep);
    n:=NormRep(resrep);
    if Length(n)=0 then
        return(rep);
    fi;
    vectors:=BaseSteinitzVectors(IdentityMat(rep.dimension,rep.field),NullspaceMat(n));
    sub:=Spin(rep,vectors.factorspace);
    if IsProjectiveRep(SubmoduleRep(rep,sub))=false then
        Error("The routine produced a basis called sub, for a submodule that is too big and not projective. It is only guaranteed to work for p-groups.");
        return;
    fi;
    return(sub);
end );

#!
#!
#!  IsIsomorphicSummand(rep1,rep2) only works when rep1 is indecomposable of dimension prime to the field characteristic. 
#!  It returns true if rep1 is isomorphic to a direct summand of rep2, false otherwise. It relies on a result of Benson
#! and Carlson. Written by Peter Webb Feb 2020.
#!
#!


InstallGlobalFunction( IsIsomorphicSummand, function(rep1,rep2)
    local hom, vecs, d;
    if rep1.dimension mod Characteristic(rep1.field) =0 then
        Error("The dimension needs to be prime to the characteristic.");
        return;
    fi;
    if rep1.dimension <> rep2.dimension then
        return false;
    fi;
    hom:=TensorProductRep(DualRep(rep1),rep2);
    vecs:=FixedQuotient(hom);
    d:=Length(vecs);
    vecs:=SafeBaseMat(Concatenation(vecs,FixedPoints@(hom)));
    return Length(vecs) > d;
end );



#
#!PrincipalIdempotent(group,prime) returns a vector in the representation space of 
#!RegularRep(group, GF(prime)) that is the block idempotent of the principal block. Thus if b is 
#!this idempotent and e denotes the vector in the regular representation corresponding to the
#!identity elements, the vector b.e is returned. Spinning it gives a basis for the principal block. 
#!This idempotent is constructed using a formula of B. Kulshammer Arch. Math 56 (1991) 313-319.
#!Code written by Peter Webb February 2020.
#


InstallGlobalFunction( PrincipalIdempotent, function(group,prime)
    local els, primepart, coeffs, pels, qels, qelspos, x, y, position; 
    els:=RightCosetReps(group, Group(Identity(group)));
    primepart:=prime^LogInt(Size(group),prime);
    coeffs:=List(els,x->0);
    pels:=ShallowCopy(els);
    for x in [1..Length(pels)] do
        if RemInt(primepart,Order(pels[x])) <>0 then
            Unbind(pels[x]);
        fi;
    od;
    #pels:=Compacted(pels);
    qels:=ShallowCopy(els);
    qelspos:=List(qels,x->1);
    for x in [1..Length(qels)] do
        if GcdInt(Order(qels[x]),prime)<>1 then
        Unbind(qels[x]);
        qelspos[x]:=0;
        fi;
    od;
    #qels:=Compacted(qels);
    for x in pels do
        for y in qels do
        position:=Position(els, x*y);
        coeffs[position]:=coeffs[position]+qelspos[position];
        od;
    od;
    coeffs:=(Length(qels)*Z(prime)^0)^-1*coeffs;
    return(coeffs);
end );

# -------------------------------- NEW FEATURES IN 1.0.0 ---------------------------------------- #


InstallGlobalFunction( IrreducibleReps, function(rep)
    local L, irreps, gens, i, rho;
    if not IsRepr(rep) then
        Error("for usage, see ?IrreducibleReps");
    fi;
    irreps := IrreducibleRepresentations( rep.group );
    gens := rep.generatorsofgroup;
    L:=[];
    for i in [1..Length( irreps )] do
        rho := GroupHomomorphismByImages(rep.group, Group( List( gens, x -> x^irreps[i] ) ) );
        Add( L, rho );
    od;
    rep.irreps := L;
    return rep;
end );


InstallGlobalFunction( IsRepr, function(rep)
    if IsRecord(rep) then
        if "isRepresentation" in Set(RecNames( rep )) then
            return rep.isRepresentation;
        else
            return false;
        fi;
    else
        return false;
    fi;
end );


InstallGlobalFunction( ReprForG, function(arg)
    # change L by images
    local G, images, rep, arrgim, rho, Ggens;
    arrgim:=arg;
    if not Length(arrgim) in [2..3] or not IsGroup(arrgim[1]) then
        Error("for usage, see ?ReprForG");
    fi;
    G := arrgim[1];
    images := arrgim[2];
    if images = [] then
        Error("Identity group encountered: creating representations is not properly set up for the identity group.");
        return;
    fi;
    if Length(arrgim) = 2 then
        rho := GroupHomomorphismByImages(G, Group( images ) );
    elif Length(arrgim) = 3 then
        Ggens := arrgim[3];
        rho := GroupHomomorphismByImages(G, Group( images ), Ggens, images );
    fi;
    return rho;
end );


InstallGlobalFunction( IrreducibleRepsOfGroup, function(G)
    local L, irreps, gens, i, rho;
    if not IsGroup(G) then
        Error("for usage, see ?IrreducibleRepsOfGroup");
    fi;
    irreps := IrreducibleRepresentations( G );
    gens := GeneratorsOfGroup(G);
    L:=[];
    for i in [1..Length( irreps )] do
        rho := GroupHomomorphismByImages( G, Group( List( gens, x -> x^irreps[i] ) ) );
        Add( L, rho );
    od;
    return L;
end );


InstallGlobalFunction( TensorProductReps, function( rep1, rep2 )

    local gens, mgens, genimages1, genimages2, rho;
    if Source( rep1 ) <> Source( rep1 ) then
        Error("You must have two representations of the same group");
        return;
    fi;
    gens := GeneratorsOfGroup( Source( rep1 ) );
    genimages1 := GeneratorsOfGroup( Images( rep1 ) );
    genimages2 := GeneratorsOfGroup( Images( rep2 ) );
    mgens := List( [1..Length( gens )],
                    i -> TensorProductMatrix( genimages1[i], genimages2[i] ) );
    rho := GroupHomomorphismByImagesNC( Source( rep1 ), Group( mgens ) );
    return rho;
end );

InstallGlobalFunction( DiagonalRep, function( rep )
    local rho;
    rho := REPN_ComputeUsingSerre(rep);
    if IsBound\.( rho, RNamObj( "diagonal_rep" ) ) then
        return rho.diagonal_rep;
    fi;
end );

# -------------------------------- END NEW FEATURES ---------------------------------------- #

GroupRepOps:=rec(Decompose:=DecomposeGroupRep,
                            SubmoduleRep:=SubmoduleGroupRep,
                            QuotientRep:=QuotientGroupRep,
                            Spin:=GroupSpin,
                            CoSpin:=GroupCoSpin,
                            HomBasis:=GroupHomBasis,
                            DecomposeSubmodule:=GroupDecomposeSubmodule,
                            SumOfImages:=GroupSumOfImages,
                            TensorProductRep:=TensorProductGroupRep,
                            DirectSumRep:=DirectSumGroupRep			
                            );