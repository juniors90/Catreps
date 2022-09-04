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
InstallGlobalFunction( EndomorphismGroup, function( cat, obj )
    local i, g, generators, permutations;    
    permutations := [ ];
    generators := StructuralCopy( cat.generators );
    for g in generators do
        if Origin( cat, g ) = obj and Terminus( cat, g ) = obj then
            # add identity
            i:=1;
            for i in [1..Length(g)] do
                if not IsBound(g[i]) then
                    g[i]:=i;
                fi;
                i:=i+1;
            od;
            # convert to permutation and add to list
            Add( permutations, PermList( g ) );
        fi;
    od;

    if Length( permutations ) > 0 then
        return Group( permutations );
    fi;

    return Group( [], () );
end );


InstallGlobalFunction( MorphismsRep, function(rep)
    local cat, n, genmat, g, mormat, oldlength, newlength, i, j, k, h, templist;
    cat := rep.category;
    if IsBound( rep.morphimgs ) then
        return( rep.morphimgs );
    fi;

    if not IsBound( cat.objects ) then
        Objects(cat);
    fi;
    
    n := Length( cat.objects );

    genmat:=EmptyMat( n, n );
    mormat:=EmptyMat( n, n );
    
    for g in [1..Length(rep.genimages)] do
        Add( genmat[ Origin( cat, cat.generators[g] ) ][ Terminus( cat, cat.generators[g]) ], rep.genimages[g] );
    od;

    for i in [1..n] do
        Add( mormat[i][i], SafeIdentityMat( rep.dimension[i], rep.field ) );
    od;
    
    oldlength := 0;
    newlength := Sum( List( mormat, x->Sum( List( x, Length ) ) ) );
    
    while oldlength < newlength do
        oldlength := newlength;
        for i in [1..n] do
            for j in [1..n] do
                templist := [];
                for k in [1..n] do
                    for g in genmat[k][j] do
                        for h in mormat[i][k] do
                            if h<>[] and g<>[] then
                                Add( templist, h * g );
                            fi;
                        od;
                    od;
                 od;
                 Append( mormat[i][j], templist );
                 mormat[i][j] := SSortedList( mormat[i][j] );
            od;
        od;
        newlength := Sum( List( mormat, x->Sum( List( x, Length ) ) ) );
    od;
    
    rep.morphimgs:=mormat;
    
    return( mormat );
end );


InstallGlobalFunction( Evaluation, function( rep, obj )
    local g, genMats;
    genMats := [];
    for g in [1..Length( rep.genimages )] do
        if Origin( rep.category, rep.category.generators[g] ) = obj and Terminus( rep.category, rep.category.generators[g] ) = obj then
            Add(genMats, rep.genimages[g]);
        fi;
    od;

    return Rep( EndomorphismGroup( rep.category, obj ), genMats, rep.field );
end );