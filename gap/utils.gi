#!
#! FixedPoints. . . . . produce a basis for the fixed points of a module
#!
#! Improvement by Dan Christensen of code originally by Peter Webb.
#!

FixedPoints@ := function(rep)
    local v, one, i, j;
	if IsBound(rep.fixedPoints) then
		return rep.fixedPoints;
	fi;
    if IsTrivial(rep.group) then
        return(IdentityMat(rep.dimension,rep.field));
    fi;
    v:=[];
    one:=One(rep.field);
    for i in [1..rep.dimension] do 
        v[i] := Concatenation(List(rep.genimages, m->m[i]));
        for j in [0..Length(rep.genimages)-1] do 
            v[i][i+rep.dimension*j] := v[i][i+rep.dimension*j] - one;
        od;
    od;
    v:=NullspaceMat(v);
    rep.fixedPoints:=v;
    return(v);
end ;
#  ---------------------------- categories ------------------

InstallGlobalFunction( Domains, function(C)
Objects(C);
C.domain:=List(C.generators,x->Origin(C,x));
C.codomain:=List(C.generators,x->Terminus(C,x));
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

InstallGlobalFunction( ExtractHom, function(vec, dimvec1, dimvec2)
    local n, hom, l, i, j, mat;
    n:=0;
    hom:=[];
    for l in [1..Length(dimvec1)] do
        mat:=[];
        for i in [1..dimvec1[l]] do
            mat[i]:=[];
            for j in [1..dimvec2[l]] do
                n:=n+1;
                mat[i][j]:=vec[n];
            od;
        od;
        Add(hom,mat);
    od;
    return(hom);
end );


# ----------------------- cateogoriesandgroups -----------------

InstallGlobalFunction( EmptyMat, function(r,s)
    local mat, i, j;
    mat := [];
    for i in [1..r] do
    mat[i] := [];
        for j in [1..s] do
            mat[i][j] := [];
        od;
    od;
    return(mat);
end );

InstallGlobalFunction( SafeIdentityMat, function(n, F)
    if n=0 then
        return [];
    else
        return IdentityMat(n, F);
    fi;
end );

# ----------------------- representations -----------------

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


InstallGlobalFunction( SafeBaseMat, function(M)
    if IsEmpty( M ) or IsEmpty( M[1] ) then
        return [];
    else
        return BaseMat(M);
    fi;
end );

InstallGlobalFunction( SpinFixedPts, function(rep,obj)
    local pt, e, replist, toSpin, i, output, bases;
    output:=[];
    e:=Evaluation(rep,obj);
    bases:=Decompose(e);
    replist:=List(bases,x->SubmoduleRep(e,x));
    for i in [1..Length(replist)] do
        if IsEmpty(Flat(replist[i].genimages)) then
            output[i]:=[];
        else
            pt:=FixedPoints@(replist[i]);
            if pt=[] then
                output[i]:=[];
            else
                toSpin:=EmptyMat(1,Length(rep.dimension))[1];
                toSpin[obj]:=pt*bases[i];
                output[i]:=Spin(rep,toSpin);
            fi;
        fi;
    od;
    return output;
end );

InstallGlobalFunction( BasisVecDims, function(n)
    local yr, i, j, bases, output;
    output:=[];
    if n>1 then
        yr:=YonedaRep(FI(n),3,GF(2));
        for i in [1..n+1] do
            bases:=SpinBasisVec(yr,i);
            output[i]:=[];
            for j in [1..Length(bases)] do
                output[i][j]:=List(bases[j],Length);
            od;
        od;
    fi;
    return output;
end );

InstallGlobalFunction( SpinBasisVec, function(rep,obj)
    local v, e, basislist, i, toSpin, output;
    output:=[];
    e:=Evaluation(rep,obj);
    basislist:=Decompose(e);
    for i in [1..Length(basislist)] do
        if IsZero(basislist) then
            output[i]:=[];
        else
            v:=basislist[i][1];
            toSpin:=EmptyMat(1,Length(rep.dimension))[1];
            toSpin[obj]:=[v];
            output[i]:=Spin(rep,toSpin);
        fi;
    od;
    return output;
end );

InstallGlobalFunction( IsDirectSum, function(summands,sum)
    local i, summandMat, bigMat;
    summandMat:=[];
    for i in [1..Length(summands)] do
        Append(summandMat,summands[i]);
    od;
    bigMat:=Concatenation(summandMat,sum);
    return RankMat(bigMat)=RankMat(summandMat) and RankMat(summandMat)=Length(summandMat) and RankMat(sum)=Length(sum) and Length(summandMat)=Length(sum);
end );

InstallGlobalFunction( TestOneFixedPt, function(rep)
    local i, j, e, replist;
    for i in [1..Length(rep.dimension)] do
        e := Evaluation(rep,i);
        replist := List( Decompose(e),x->SubmoduleRep(e,x));
        for j in [1..Length(replist)] do
            if not IsEmpty(Flat(replist[j].genimages)) and Length( FixedPoints@(replist[j]))>1 then
                return false;
            fi;
        od;
    od;
    return true;
end );

InstallGlobalFunction( SafeFixedPoints, function(rep)
    if rep.dimension=0 then
        return [];
    else
        return FixedPoints@(rep);
    fi;
end );

InstallGlobalFunction( BrauerTraceImage, function( rep, p )
    local resp, resh, maxsub, image, h, M;
    resp := RestrictedRep( rep.group, p, rep );
    maxsub := MaximalSubgroups(p);
    image := [];
    for h in maxsub do
        M := RelativeTrace( p, h, resp );
        resh := RestrictedRep( resp.group, h, resp );
        Append( image, SafeMatrixMult( FixedPoints@(resh), M, rep.dimension ) );
    od;
    return(SafeBaseMat(image));
end );

InstallGlobalFunction( BrauerRep, function(rep,p)
    local n, resn, resp, fp, fprep, maxsub, traceimages, h, M, resh, v, b;
    n := Normalizer( rep.group, p );
    resn := RestrictedRep( rep.group, n, rep );
    resp := RestrictedRep( rep.group, p, rep );
    fp := FixedPoints@( resp );
    fprep := SubmoduleRep( resn, fp );
    maxsub := MaximalSubgroups( p );
    traceimages := [];
    for h in maxsub do
        M:=RelativeTrace(p,h,resp);
        resh:=RestrictedRep(resp.group,h,resp);
        Append(traceimages,SafeMatrixMult(FixedPoints@(resh),M,rep.dimension));
    od;
    traceimages := SafeBaseMat( traceimages );
    v := VectorSpace( rep.field, fp );
    b := Basis(v,fp);
    traceimages := List( traceimages, x -> Coefficients( b, x ) );
    return QuotientRep( fprep, traceimages );
end );

InstallGlobalFunction( AllSpinDims, function( n, evalObj, outputRec )
    local yr, dim, i, v, output, spunDims, toSpin;
    yr := YonedaRep( FI( n ), 3, GF( 2 ) );
    dim := Evaluation( yr, evalObj + 1 ).dimension;
    if outputRec then
        output := rec();
        for v in GF(2)^dim do
            toSpin := EmptyMat(1, n+1)[1];
            toSpin[evalObj+1] := [v];
            spunDims := String( List ( Spin (yr, toSpin ), Length ) );
            if IsBound( output.(spunDims) ) then
                Append( output.(spunDims), [v] );
            else
                output.(spunDims) := [v];
            fi;
        od;
        return output;
    else
        output := [];
        i := 0;
        for v in GF(2)^dim do
            i := i + 1;
            toSpin := EmptyMat(1,n+1)[1];
            toSpin[evalObj + 1] := [v];
            output[i] := List( Spin( yr, toSpin ), Length );
        od;
        return SSortedList( output );
    fi;
end );

InstallGlobalFunction( BasisVecBases, function( n )
    local yr, i, j, bases, output;
    output := [];
    if n > 1 then
        yr := YonedaRep( FI( n ), 3, GF(2) );
        for i in [1..n+1] do
            output[i] := SpinBasisVec( yr, i );
        od;
    fi;
    return output;
end );

InstallGlobalFunction( FixedPtDims, function( n )
    local yr, i, j, bases, output;
    output := [];
    if n> 1 then
        yr := YonedaRep(FI(n),3,GF(2));
        for i in [1..n+1] do
            bases := SpinFixedPts(yr,i);
            output[i] := [];
            for j in [1..Length( bases )] do
                if bases[j] = [] then
                    output[i][j] :=[];
                else
                    output[i][j] := List( bases[j], Length );
                fi;
            od;
        od;
    fi;
    return output;
end );

InstallGlobalFunction( DimSummands, function(n, obj)
    local summands, output, i, j, k;
    summands := FISummandEvalReps(n,obj,GF(2));
    output := [];
    for i in [1..Length(summands)] do
        output[i] := [];
        for j in [1..Length(summands[i])] do
            output[i][j] := [];
            for k in [1..Length(summands[i][j])] do
                output[i][j][k] := summands[i][j][k].dimension;
            od;
        od;
    od;
    return output;
end );

InstallGlobalFunction( SafeDimHom, function( d, rep )
    if rep.dimension = 0 then
        return 0;
    else
        return DimHom( d, rep );
    fi;
end );

InstallGlobalFunction( ProjSummands, function( n, obj )
    local summands, output, i, j, k;
    summands:=FISummandEvalReps(n,obj,GF(2));
    output:=[];
    for i in [1..Length(summands)] do
        output[i]:=[];
        for j in [1..Length(summands[i])] do
            output[i][j]:=[];
            for k in [1..Length(summands[i][j])] do
                if IsEmpty(Flat(summands[i][j][k].genimages)) then
                    output[i][j][k]:=fail;
                else
                    output[i][j][k]:=IsProjectiveRep(summands[i][j][k]);
                fi;
            od;
        od;
    od;
    return output;
end );

InstallGlobalFunction( SocleSeries, function(rep)
    local radseries, n, temp, socseries, i;
    radseries:=RadicalSeries(DualRep(rep));
    socseries:=[];
    socseries[1]:=[];
    socseries[2]:=[];
    temp:=[];
    temp[1]:=List(radseries[1],y->SocleNullspaceMat(TransposedMat(y),rep.dimension,rep.field));
    temp[2]:=List(radseries[2],DualRep);
    n:=Length(radseries[2]);
    for i in [1..n+1] do
        socseries[1][i]:=temp[1][n+2-i];
    od;
    for i in [1..n] do
        socseries[2][i]:=temp[2][n+1-i];
    od;
    return(socseries);
end );

InstallGlobalFunction( RemoveFromBottom, function( rep, d )
    local newrep;
    newrep := QuotientRep( rep, SumOfImages( d, rep ) );
    return newrep;
end );

InstallGlobalFunction( RadicalSeries, function(rep)
    local rad, i, submod, submodrad, layer;
    rad:=[];
    layer:=[];
    rad[1]:=SafeIdentityMat(rep.dimension,rep.field);
    i:=1;
    while Length(rad[i])>0 do
        i:=i+1;
        submod:=SubmoduleRep(rep,rad[i-1]);
        submodrad:=RadicalRep(submod);
        rad[i]:=SafeMatrixMult(submodrad,rad[i-1],submod.dimension);
        layer[i-1]:=QuotientRep(submod,submodrad);
    od;
    return([rad,layer]);
end );

InstallGlobalFunction( DisplayButterflyDims, function(sn, subGens)
    local m, b;
    m:=PermutationRepOnCosets(sn,Subgroup(sn,subGens),GF(2));
    b:=ButterflyFactors(m,SocleSeries(m)[1],RadicalSeries(m)[1]);
    DisplayMatrix(List(b,x->List(x,y->y.dimension)));
    return b;
end );

InstallGlobalFunction( ButterflyDimsRep, function(rep)
    local b;
    b := ButterflyFactors( rep, SocleSeries( rep )[1], RadicalSeries( rep )[1] );
    DisplayMatrix(List(b,x->List(x,y->y.dimension)));
    return b;
end );

InstallGlobalFunction( npButterflyDimsRep, function(rep)
    return ButterflyFactors(rep,SocleSeries(rep)[1],RadicalSeries(rep)[1]);
end );

InstallGlobalFunction( DisplayMatrix, function( A )
    PrintArray( List (A, x -> List( x, Int ) ) );
end );

InstallGlobalFunction( ExamineButterflyFactors, function(b,dRec)
    local name;
    Print("Original matrix: \n");
    DisplayMatrix(List(b,x->List(x,y->y.dimension)));
    Print("\n\nFixed points:\n");
    DisplayMatrix(List(b,x->List(x,y->Length(SafeFixedPoints(y)))));
    Print("\n\nFactors:\n");
    for name in RecNames(dRec) do
        Print(name);
        Print("\n");
        DisplayMatrix(List(b,x->List(x,y->SafeDimHom(dRec.(name),y))));
        Print("\n");
    od;
end );

InstallGlobalFunction( FISummandEvalReps, function(n,obj)
    local catreplist, yr, i, e, j, output;
    output:=[];
    yr:=YonedaRep(FI(n),obj+1,GF(2));
    catreplist:=List(Decompose(yr),x-> SubmoduleRep(yr,x));
    for i in [1..Length(catreplist)] do
        output[i]:=[];
        for j in [1..n+1] do
            e:=Evaluation(catreplist[i],j);
            output[i][j]:=List(Decompose(e),x-> SubmoduleRep(e,x));
        od;
    od;
    return output;
end );

# ------------------------- vectors ------------------------------------- #
InstallGlobalFunction( OrthogonalComplement, function( veclist, dim, field )
    if veclist = [] then
        return( IdentityMat( dim, field ) );
    else
        return NullspaceMat( TransposedMat( veclist ) );
    fi;
end );

InstallGlobalFunction( GeneratorDomains, function(cat)
    local M, i;
    if IsBound(cat.generatordomains) then
        return cat.generatordomains;
    fi;
    M:=EmptyMat(Length(cat.objects),Length(cat.objects));
    for i in [1..Length(cat.generators)] do
        Add(M[cat.domain[i]][cat.codomain[i]],i);
    od;
    cat.generatordomains:=M;
    return M;
end );

# -------------------------- homomorphisms ------------------

InstallGlobalFunction( RemoveFromTop, function ( rep, d )
    local  newrep;
    newrep := SubmoduleRep( rep, KernelIntersection( rep, d ) );
    return newrep;
end );

InstallGlobalFunction( KernelIntersection, function( g, h )
    local transposedmats;
    transposedmats := List( HomBasis( g, h ), TransposedMat);
    if Length(transposedmats) = 0 then 
        return( SafeIdentityMat( g.dimension, g.field ) );
    else
        return SafeNullspaceMat( TransposedMat( Concatenation(transposedmats ) ) , g.field );
    fi;
end );

InstallGlobalFunction( RadicalRep, function(rep)
    if rep.dimension=0 then return([]); fi;
    return(MTX.BasisRadical(RepToMeataxeModule(rep)));
end );


InstallGlobalFunction( RepToMeataxeModule, function(rep)
    if rep.genimages=[] then
        return(GModuleByMats([], rep.dimension, rep.field ));
    elif rep.dimension=0 then
        return(GModuleByMats([], rep.dimension, rep.field ));
    else
        return(GModuleByMats(rep.genimages, rep.field ));
    fi;
end );




InstallGlobalFunction( DualRep, function(rep)
    if rep.dimension=0 then
        return(rep);
    fi;
    return rec(
     group:=rep.group,
     genimages:=List(InverseGenImages(rep), TransposedMat),
     field:=rep.field,
     dimension:=rep.dimension,
     isRepresentation:=true,
     operations:=CatRepOps
     );
end );

InstallGlobalFunction( RestrictedRep, function(G,H,rep)
    local R;
    R:=MatricesOfElements(rep,GeneratorsOfGroup(H));
    return rec(
        group:=H,
        genimages:=R,
        field:=rep.field,
        dimension:=rep.dimension,
        isRepresentation:=true,
        operations:=CatRepOps
        );
end );

InstallGlobalFunction( MatricesOfElements, function(rep,l)
local newrep, newl, matlist, s, temp, x, n;
    MakeIsoToPermGroup(rep);
    newrep:=rec(
     group:=rep.permgroup,
     genimages:=rep.genimages,
     field:=rep.field,
     dimension:=rep.dimension,
     isRepresentation:=true,
     operations:=CatRepOps
     );
    n:=Length(l);
    newl:=List(l, x->Image(rep.isotopermgroup,x));
    matlist:=List(l,x->IdentityMat(rep.dimension,rep.field));
    while Size(newrep.group)>1 do
        s:=SmallestMovedPoint(newrep.group);
        temp:=MatsOfCosetReps(newrep);
        newrep:=temp[1];
        matlist:=List([1..n], x->temp[2][s^newl[x]]*matlist[x]);
        newl:=List(newl,x->x*temp[3][s^x]^-1);
    od;
    return(matlist);
end );

InstallGlobalFunction( RelativeTrace, function(P,Q,rep)
   local M,l,g,ghbi;
   M:=NullMat(rep.dimension,rep.dimension,rep.field);
   l:=RightCosetReps(P,Q);
   ghbi:=RepToGHBI(rep );
   for g in l do
      M:=M+Image(ghbi,g);
   od;
   return M;
end );




InstallGlobalFunction( IsProjectiveRep, function(rep)
    local s, resrep, n;
    if Characteristic(rep.field)=0 then return(true);
    fi;
    s:=SylowSubgroup(rep.group, Characteristic(rep.field));
#    s:=Subgroup(rep.group,SmallGeneratingSet(s));  This line is for permutation groups
    resrep:=RestrictedRep(rep.group, s, rep);
    n:=NormRep(resrep);
    if Rank(n)*Size(s) = rep.dimension then
        return(true);
    else return(false);
    fi;
end );

InstallGlobalFunction( SocleNullspaceMat, function(M,n,F)
    if n=0 or IsEmpty(M) then
        return SafeIdentityMat(n,F);
    else
        return NullspaceMat(M);
    fi;
end );

InstallGlobalFunction( InverseGenImages, function(rep)
    if IsBound(rep.inversegenimages) then
        return(rep.inversegenimages);
    fi;    
    rep.inversegenimages:=List(rep.genimages, x->x^-1);
    return(rep.inversegenimages);
end );

InstallGlobalFunction( OldInverseGenImages, function(rep)
    local orders, inverses, i;
    if IsBound(rep.inversegenimages) then
        return(rep.inversegenimages);
    fi;    
    orders:=List(GeneratorsOfGroup(rep.group),Order);
    inverses:=[];
    for i in [1..Length(orders)] do
        if orders[i]=infinity then
            inverses[i]:=rep.genimages[i]^-1;
        else
            inverses[i]:=rep.genimages[i]^(orders[i]-1);
        fi;
    od;
    rep.inversegenimages:=inverses;
    return(inverses);
end );

InstallGlobalFunction( MakeIsoToPermGroup, function(rep)
    if IsBound(rep.permgroup) then
        return;
    fi;
    rep.isotopermgroup:=IsomorphismPermGroup(rep.group);
    rep.permgroup:=Image(rep.isotopermgroup,rep.group);
    return;
end );

InstallGlobalFunction( MatsOfCosetReps, function(rep)
    local s, orbit, positions, mats, inversemats, perms, gens, q, x, j, ximage, table,
    subgroupgens, subgroupgenimages;
    s:=SmallestMovedPoint(rep.group);
    orbit:=[s];
    positions:=[];
    mats:=[];
    mats[s]:=IdentityMat(rep.dimension,rep.field);
    inversemats:=[];
    inversemats[s]:=IdentityMat(rep.dimension,rep.field);
    perms:=[];
    perms[s]:=();
    gens:=GeneratorsOfGroup(rep.group);
    InverseGenImages(rep);
    table:=[];
    for j in [1..Length(gens)] do
        table[j]:=[];
    od;
    q:=1;
    while IsBound(orbit[q]) do
        x:=orbit[q];
        for j in [1..Length(gens)] do
            ximage:=x^gens[j];
            if not ximage in orbit then
                Add(orbit,ximage);
                positions[ximage]:=j;
                perms[ximage]:=perms[x]*gens[j];
                mats[ximage]:=mats[x]*rep.genimages[j];
                inversemats[ximage]:=rep.inversegenimages[j]*inversemats[x];
                else table[j][x]:=perms[x]*gens[j]*perms[ximage]^-1;
            fi;
        od;
        q:=q+1;
    od;
    subgroupgens:=[];
    subgroupgenimages:=[];
    for j in [1..Length(gens)] do
        for x in orbit do
            if IsBound(table[j][x]) and
            Size(Group(subgroupgens,()))<
            Size(Group(Concatenation(subgroupgens,[table[j][x]])))
                then Add(subgroupgens,table[j][x]);
                Add(subgroupgenimages,mats[x]*rep.genimages[j]*inversemats[x^gens[j]]);
            fi;
        od;
    od;
    return ([rec(
        group:=Group(subgroupgens,()),
        genimages:=subgroupgenimages,
        field:=rep.field,
        dimension:=rep.dimension,
        isRepresentation:=true,
        operations:=CatRepOps
        ),
        mats,perms]);
end );

InstallGlobalFunction( RightCosetReps, function(G,H)
   return List(RightCosets(G,H), Representative);
end );


InstallGlobalFunction( RepToGHBI, function(phi)
    local glg,h;
    glg := GLG(phi.dimension,Size(phi.field));
    h := GroupHomomorphismByImagesNC(
            phi.group,glg,
            GeneratorsOfGroup(phi.group),
            phi.genimages
        );
    #   h.isMapping:=true;
    #   h.isHomomorphism:=true;
    #   h.isGroupHomomorphism:=true;
    return h; 
end );


InstallGlobalFunction( NormRep, function(rep)
    local output, newrep, temp, a, sum;
    output:=IdentityMat(rep.dimension,rep.field);
    MakeIsoToPermGroup(rep);
    newrep:=rec(
        group:=rep.permgroup,
        genimages:=rep.genimages,
        field:=rep.field,
        dimension:=rep.dimension,
        isRepresentation:=true,
        operations:=CatRepOps
    );
    while Size(newrep.group)>1 do
        temp:=MatsOfCosetReps(newrep);
        newrep:=temp[1];
        sum:=NullMat(rep.dimension,rep.dimension,rep.field);
        for a in temp[2] do
            sum:=sum+a;
        od;
        output:=sum*output;
    od;
    return output;
end );


InstallGlobalFunction( OldNormRep, function(rep)
        local p, currentgroup, stablist, subgroup, replist, norm;
        currentgroup := rep.group;
        p := RepToGHBI(rep);
        norm := IdentityMat(rep.dimension,rep.field);
        while not IsTrivial(currentgroup) do;
             stablist := StabChain(currentgroup);
             subgroup := Subgroup(currentgroup,stablist.stabilizer.generators);
             replist := RightCosetReps(currentgroup, subgroup);
             replist := List( replist, x -> ImagesRepresentative(p, x));
             norm := Sum(replist) * norm;
             currentgroup := subgroup;
        od;
        return norm;
end );

InstallGlobalFunction( GLG, function(n,q)
   if n=1 then
      return(Group(Z(q)*IdentityMat(1,GF(q))));
   else 
      return(GeneralLinearGroup(n,q));
   fi;
end );


InstallGlobalFunction( TensorProductRep, function(g,h)
    return g.operations.TensorProductRep( g, h );
end );



CatRepOps:=rec(Decompose:=DecomposeCatRep,
                            SubmoduleRep:=SubmoduleCatRep,
                            QuotientRep:=QuotientCatRep,
                            Spin:=CatSpin,
                            CoSpin:=CatCoSpin,
                            HomBasis:=CatHomBasis,
                            DecomposeSubmodule:=CatDecomposeSubmodule,
                            SumOfImages:=CatSumOfImages,
                            TensorProductRep:=TensorProductCatRep,
                            DirectSumRep:=DirectSumCatRep
                            );