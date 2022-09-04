# ----------------------------------------------------------------------

InstallGlobalFunction( Spin, function(rep,veclist)
    return rep.operations.Spin(rep,veclist);
end );

InstallGlobalFunction( CatSpin, function(rep,veclistlist)
    local basis, olddim, dim, newlist, g, n, v;
    basis:=List(veclistlist,SafeBaseMat);
    dim:=Sum(List(basis,Length));
    olddim:=dim-1;
    while dim>olddim do
        olddim:=dim;
        newlist:=List(basis,x->List(x,y->y));
        for n in [1..Length(rep.genimages)] do
            g:=rep.genimages[n];
            for v in basis[rep.category.domain[n]] do
                Add(newlist[rep.category.codomain[n]],v*g);
            od;
        od;
        basis:=List(newlist,SafeBaseMat);
        dim:=Sum(List(basis,Length));
    od;
    return basis;
end );

#! CoSpin(rep, list of lists of vectors)
#! returns a basis for the largest subrepresentation contained in the
#! subspaces spanned by the vectors in the lists.
#! There is one entry in the list (of lists of vectors) for each object in the category,
#! and it is a list of vectors which all belong to the representation space at that object.
#!
#! CatCoSpin is called by CoSpin(rep, list of lists of vectors) when rep is a
#! category representation.
#!
InstallGlobalFunction( CoSpin, function(rep, veclist)
    return rep.operations.CoSpin(rep, veclist);
end );

InstallGlobalFunction( CatCoSpin, function( rep, veclistlist )
    local basis, olddim, dim, newlist, g, n, v, output, i, transposedimages;
    transposedimages := List( rep.genimages, TransposedMat );
    basis := [];
    for i in [1..Length(rep.dimension)] do
        basis[i] := OrthogonalComplement( veclistlist[i], rep.dimension[i], rep.field );
    od;
    dim := Sum(List(basis,Length));
    olddim := dim-1;
    while dim > olddim do
        olddim := dim;
        newlist := List(basis,x->List(x,y->y));
        for n in [1..Length(rep.genimages)] do
            g := transposedimages[n];
            for v in basis[rep.category.codomain[n]] do
                Add( newlist[ rep.category.domain[n] ], v * g );
            od;
        od;
        basis := List( newlist, SafeBaseMat );
        dim := Sum( List( basis, Length ) );
    od;
    output := [];
    for i in [1..Length( rep.dimension )] do
        output[i] := OrthogonalComplement( basis[i], rep.dimension[i], rep.field );
    od;
    return output;
end );

#! SumOfImages(rep,rep) . . returns a basis for the sum of images of all 
#! module homomorphisms A -> B
#!
#! GroupSumOfImages is called by SumOfImages(rep,rep) when rep is a
#! group representation.
#!
#! Updated by Peter Webb Sep 11, 2007
#!

InstallGlobalFunction( SumOfImages, function(M,N)
    return M.operations.SumOfImages(M,N);
end );

InstallGlobalFunction( CatSumOfImages, function(M,N)
    local  homs, basis, l, h;
    homs:=HomBasis(M,N);
    basis:=[];
    for l in [1..Length(M.dimension)] do
        basis[l]:=[];
        for h in homs do
            Append(basis[l],h[l]);
        od;
        if not(basis[l]=[]) then
            basis[l]:=List(SafeBaseMat(basis[l]),x->x);
        fi;
    od;
    return(basis);
end );

InstallGlobalFunction( Decompose, function(rep)
    return rep.operations.Decompose(rep);
end );


InstallGlobalFunction( DecomposeCatRep, function(rep)
	local summands, result, q;
	if IsBound(rep.summands) then return(rep.summands);
	fi;
	summands := [List(rep.dimension, x->SafeIdentityMat(x,rep.field))];
	q := 1;
	# We maintain the following invariants:
	# - summands is a list of basis structures; the 'union' of these
	#   lists forms a basis structure for rep.
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

InstallGlobalFunction( DecomposeSubmodule, function(M, basisstructure)
    return M.operations.DecomposeSubmodule(M, basisstructure);
end );

InstallGlobalFunction( CatDecomposeSubmodule, function(M,basisstructure)
	local newrep, c, d, initlist, b, x, a, n, z, kernel, kdim, image ;
	newrep := SubmoduleRep(M, basisstructure);
	d:=Maximum(newrep.dimension);
	initlist := HomBasis(newrep,newrep);
	Add(initlist, List(initlist[1], u->u*0));
	b := Length(initlist);
	while b > 0 do;
		c := b - 1;
		while c > 0 do;
	        a:=initlist[b]+initlist[c];
			n:=1;
			while n < d do
				for z in [1..Length(a)] do
		    		if not(a[z]=[]) then
		        		a[z]:=a[z]*a[z];
		    		fi;
				od;
				n:=2*n;
			od;
			kernel:=List(a,u->SafeNullspaceMat(u,M.field));
		 	kdim:=Sum(List(kernel,Length));
			if not(kdim=0 or kdim=Sum(newrep.dimension)) then
				image:=List(a,SafeBaseMat);
				for z in [1..Length(a)] do
				    if image[z] <> [] then
				        image[z]:=image[z]*basisstructure[z];
				    fi;
				    if kernel[z]<>[] then
				        kernel[z]:=kernel[z]*basisstructure[z];
				    fi;
				od;
				return [kernel,image];
			fi;
	   		c := c - 1;
		od;
		b := b - 1;
	od;
	return [basisstructure];
end );


InstallGlobalFunction( SubConstant, function(rep)
    local morphisms, vecListList, obj1, obj2, m, m1, m2, perpMat, perpSum, vec;
    morphisms:=MorphismsRep(rep);
    vecListList:=[];
    if Length(morphisms)>0 then
        for obj1 in [1..Length(morphisms)] do
            perpSum:=[];

            for obj2 in [1..Length(morphisms)] do
                perpMat:=[];

                if obj1<>obj2 and Length(morphisms[obj1][obj2])>1 then
                    for m1 in morphisms[obj1][obj2] do
                        for m2 in morphisms[obj1][obj2] do
                            if m1<>m2 then
                                for vec in TransposedMat(m1-m2) do
                                    Add(perpMat,vec);
                                od;
                            fi;
                        od;
                    od;
                    Append(perpSum, SafeBaseMat(perpMat));

                elif obj1=obj2 then
                    for m in morphisms[obj1][obj1] do
                        if rep.dimension[obj1]<>0 then
                            for vec in TransposedMat(m-SafeIdentityMat(rep.dimension[obj1],rep.field)) do
                                Add(perpMat,vec);
                            od;
                        fi;
                    od;
                    Append(perpSum, SafeBaseMat(perpMat));
                fi;
            od;
            
            if perpSum<>[] then
            vecListList[obj1]:=SafeNullspaceMat(TransposedMat(perpSum),rep.field);
            else vecListList[obj1]:=SafeIdentityMat(rep.dimension[obj1],rep.field);
            fi;
        od;
    fi;

    return(vecListList);
end );