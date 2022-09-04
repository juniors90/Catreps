
InstallGlobalFunction( CatRep, function( C, L, F)
    local dimvec, i;
    dimvec := List( Objects( C ), x -> 0 );

    for i in [1..Length( C.generators )] do
        if IsBound( L[i] ) then
            dimvec[C.domain[i]] := Length( L[i] );
        fi;
    od;
    return rec(
        category := C,
        genimages := L,
        field := F,
        dimension := dimvec,
        operations := CatRepOps
    );
end );


InstallGlobalFunction( YonedaRep, function( cat, object, F)
    local genmatrices, mor, dimvec, i, j, k, f, matrix;
    genmatrices := [];
    mor := Morphisms( cat );
    dimvec := List( mor[object], Length );
    for i in [1..Length(cat.generators)] do
        if dimvec[cat.domain[i]] > 0 then
            matrix :=List( [1..dimvec[cat.domain[i]]], x -> List( [1..dimvec[cat.codomain[i]]] , y -> Zero( F ) ) );
            for j in [1..Length(mor[object][cat.domain[i]])] do
                f := Composition( mor[object][cat.domain[i]][j], cat.generators[i]);
                k := Position( mor[object][cat.codomain[i]], f );
                matrix[j][k] := One(F);
            od;
            genmatrices[i] := matrix;
        else
            genmatrices[i] := [];
        fi;
    od;
    return rec(
        category := cat,
        genimages := genmatrices,
        field := F,
        dimension := dimvec,
        operations := CatRepOps
    );
end );

#! YonedaDualRep(category, object number, field)
#! This returns the representation of the category on the dual of the subspace of the
#! category algebra spanned by the morphisms whose codomain is the specified
#! object. The representation is injective, because its dual is projective, by Yoneda's lemma.
InstallGlobalFunction( YonedaDualRep, function( cat, object, F )
    local genmatrices, mor, dimvec, i, j, k, f, matrix;
    genmatrices := [];
    mor := Morphisms(cat);
    dimvec := List( mor, x->Length( x[object] ) );
    for i in [1..Length( cat.generators )] do
        if dimvec[ cat.domain[i] ] > 0 then
            matrix := List( [1..dimvec[cat.domain[i]]], x -> List( [1..dimvec[cat.codomain[i]]], y -> Zero( F ) ) );
            for j in [1..dimvec[cat.codomain[i]]] do
                f := Composition( cat.generators[i], mor[cat.codomain[i]][object][j] );
                k := Position( mor[cat.domain[i]][object], f );
                matrix[k][j] := One(F);
            od;
            genmatrices[i] := matrix;
        else genmatrices[i] := [];
        fi;
    od;
    return rec(
        category := cat,
        genimages := genmatrices,
        field := F,
        dimension := dimvec,
        operations := CatRepOps
    );
end );


InstallGlobalFunction( DirectSumCatRep, function( rep1, rep2 )
    local i, row, col, newmat, gimages, dDim1, dDim2, cDim1, cDim2, cat;
    if rep1.field <> rep2.field or rep1.category <> rep2.category then
        Error("Representations incompatible.");
    fi;
    cat := rep1.category;
    gimages := [];
    for i in [1..Length(cat.generators)] do
        dDim1 := rep1.dimension[cat.domain[i]];
        dDim2 := rep2.dimension[cat.domain[i]];
        cDim1 := rep1.dimension[cat.codomain[i]];
        cDim2 := rep2.dimension[cat.codomain[i]];
        newmat := EmptyMat(dDim1+dDim2,cDim1+cDim2);
        if newmat <> [] and newmat[1] <> [] then
            for row in [1..Length(newmat)] do
                for col in [1..Length(newmat[1])] do
                    if row in [1..dDim1] and col in [1..cDim1] then
                        newmat[row][col] := rep1.genimages[i][row][col];
                    elif row in [dDim1+1..dDim1+dDim2] and col in [cDim1+1..cDim1+cDim2] then
                        newmat[row][col] := rep2.genimages[i][row-dDim1][col-cDim1];
                    else
                        newmat[row][col] := Zero(rep1.field);
                    fi;
                od;
            od;
        fi;
        Add( gimages, newmat );
    od;
    return rec(
        category:=cat,
        genimages:=gimages,
        field:=rep1.field,
        dimension:=rep1.dimension+rep2.dimension,
        operations:=CatRepOps
    );
end );


InstallGlobalFunction( TensorProductCatRep, function(g,h)
    local mgens, i;
    mgens:=[];
    for i in [1..Length(g.genimages)] do
        if IsBound(g.genimages[i]) and IsBound(h.genimages[i]) then
            mgens[i]:=TensorProductMatrix(g.genimages[i],h.genimages[i]);
        fi;
    od;
    return CatRep(g.category, mgens, g.field);
end);

#! SubmoduleGroupRep is called by SubmoduleRep(rep, list of vecs) when rep is a
#! group representation.

InstallGlobalFunction( SubmoduleRep, function(rep,v)
    return rep.operations.SubmoduleRep( rep, v );
end );

InstallGlobalFunction( SubmoduleCatRep, function( rep, v )
    local dimvec, vs, base, i, newimages, mat, domain, codomain;
    domain := rep.category.domain;
    codomain := rep.category.codomain;
    dimvec := List(v,Length);
    vs := [];
    base := [];
    for i in [1..Length(v)] do
        if not(v[i]=[]) then
            vs[i] := VectorSpace( rep.field, v[i] );
            base[i] := Basis( vs[i],v[i] );
        fi;
    od;
    newimages := [];
    for i in [1..Length( rep.genimages )] do
        mat := List( v[domain[i]], x-> []);
        if not( v[codomain[i]] = [] or v[domain[i]] = [] ) then
            mat := List( base[domain[i]], 
                         b -> Coefficients( base[codomain[i]], b * rep.genimages[i] ) );
        fi;
        Add(newimages, mat);
    od;
    return rec(
        category := rep.category,
        genimages := newimages,
        field := rep.field,
        dimension := dimvec,
        isRepresentation := true,
        operations := CatRepOps
    );
end );

#! QuotientGroupRep is called by QuotientRep(rep, list of vecs) when rep is a
#! group representation.

InstallGlobalFunction( QuotientRep, function(rep,v)
    return rep.operations.QuotientRep( rep, v);
end );

InstallGlobalFunction( QuotientCatRep, function(rep,v)
    local basestructure, d, zero,onemat, i, j, positions, b, p, g,
          mat, newb, newimages, baseinverse;
    if Length( v ) = 0 then
        return rep;
    fi;
    basestructure := List( v, x -> ShallowCopy( SafeBaseMat( x ) ) );
    for g in basestructure do
        TriangulizeMat(g);
    od;
    d := List( basestructure, Length);
    if d = rep.dimension then
        return ZeroCatRep( rep.category, rep.field );
    fi;
    zero := Zero(rep.field);
    positions := [];
    baseinverse := [];
    for j in [1..Length( rep.dimension )] do
        onemat:=IdentityMat( rep.dimension[j], rep.field );
        i:=1;
        positions[j] := [];
        for b in basestructure[j] do
            while b[i] = zero do
                Add( positions[j], i );
                i := i+1;
            od;
            i := i + 1;
        od;
        Append( positions[j], [i..rep.dimension[j]] );
        for p in positions[j] do
            Add( basestructure[j], onemat[p] );
        od;
        if basestructure[j] = [] then
            baseinverse[j] := [];
        else
            baseinverse[j] := basestructure[j]^-1;
        fi;
    od;
    newimages := [];
    for g in [1..Length( rep.genimages )] do
        mat:=[];
        for p in positions[rep.category.domain[g]] do
            if baseinverse[rep.category.codomain[g]] <> [] then
                b:=rep.genimages[g][p]*baseinverse[rep.category.codomain[g]];
            fi;
            newb:=[];
            for i in [d[rep.category.codomain[g]]+1..rep.dimension[rep.category.codomain[g]]] do
                Add( newb, b[i] );
            od;
            Add( mat, newb );
        od;
        Add( newimages, mat );
    od;
    return rec(
        category := rep.category,
        genimages := newimages,
        field := rep.field,
        dimension := rep.dimension-d,
        isRepresentation := true,
        operations := CatRepOps
    );
end );


#! ZeroCatRep(cat,field) returns the zero representation of the category.
InstallGlobalFunction( ZeroCatRep, function( cat, F )
    local dimvec;
    dimvec := List( Objects( cat ), x -> 0 );
    return rec(
        category := cat,
        genimages := List( cat.generators, x->[] ),
        field := F,
        dimension := dimvec,
        operations := CatRepOps
    );
end );

#! ConstantRep(cat,field) returns the constant (or trivial)
#! representation of the category.
InstallGlobalFunction( ConstantRep, function(cat,F)
    local dimvec;
    dimvec := List(Objects(cat),x->1);
    return rec(
        category := cat,
        genimages := List( cat.generators, x->[[One(F)]] ),
        field := F,
        dimension := dimvec,
        operations := CatRepOps);
end );

# ------------------------------------

InstallGlobalFunction( FixedPtBases, function(n)
    local yr, i, j, bases, output;
    output:=[];
    if n>1 then
        yr := YonedaRep( FI(n), 3, GF(2));
        for i in [1..n+1] do
            bases := SpinFixedPts( yr, i );
            output[i] := [];
            for j in [1..Length(bases)] do
                if not bases[j]=[] then
                    output[i] := bases[j];
                    break;
                fi;
            od;
        od;
    fi;
    return output;
end );

#! TestYRSummand(basis,eval,summandi) takes in a basis, an evaluation of a Yoneda Rep over
#! GF(2), and the index of the summand of that evaluation that the basis is thought to be
#! equivalent to; it returns whether it is in fact equivalent.
#!
#! Written by Moriah Elkin Spring 2019.
InstallGlobalFunction( TestYRSummand, function(basis,eval,summandi)
    local summands;
    summands := Decompose(eval);
    summands[summandi] := basis;
    return IsDirectSum( summands, BasisVectors( Basis( GF(2)^Length( basis[1] ) ) ) );
end );

