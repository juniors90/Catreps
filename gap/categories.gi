InstallGlobalFunction( Rep, function(arg)
    local G, L, rep;
    G:=arg[1];L:=arg[2];
    if L=[] then 
        Error("Identity group encountered: creating representations is not properly set up for the identity group.");
        return;    
    fi;
    rep:=rec(
            group:=G,
            genimages:=L,
            isRepresentation:=true,
            dimension:=Length(L[1]),
            operations:=CatRepOps
            );
    if Length(arg)=2 then
        rep.field:=Field(L[1][1]);
    else
        rep.field:=arg[3];
    fi;
    return rep;
end );

InstallGlobalFunction( SupportOfMorphism, function(m)
    local s, i;
    s:=[];
    for i in [1..Length(m)] do
        if IsBound(m[i]) then Add(s,i);
        fi;
    od;
    return(s);
end );

InstallGlobalFunction( IdentityMorphism, function(l)
    local m, i;
    m:=[];
    for i in l do
        m[i]:=i;
    od;
    return(m);
end );

InstallGlobalFunction( Composition, function(f,g)
    local i, h;
    h:=[];
    for i in [1..Length(f)] do
        if IsBound(f[i]) then
            if IsBound(g[f[i]]) then
                h[i]:=g[f[i]];
            else return(false);
            fi;
        fi;
    od;
    return(h);
end );

InstallGlobalFunction( IsComposable, function(f,g)
    local i;
    for i in [1..Length(f)] do
        if IsBound(f[i]) then
            if not IsBound(g[f[i]]) then return(false);
            fi;
        fi;
    od;
    return(true);
end );

InstallGlobalFunction( Objects, function( cat )
    local m;
    if IsBound( cat.objects ) then
        return(cat.objects);
    fi;
    cat.objects := [];
    for m in cat.generators do
        Add( cat.objects, SupportOfMorphism(m) );
    od;
    cat.objects := SSortedList( cat.objects );
    return( cat.objects );
end );

InstallGlobalFunction( Origin, function(cat,m)
    local k, x;
    k:=PositionBound(m);
    for x in [1..Length(cat.objects)] do
         if k in cat.objects[x] then
              return(x);
         fi;
    od;
end );

InstallGlobalFunction( Terminus, function( cat, m )
    local k, x;
    k := m[PositionBound(m)];
    for x in [1..Length(cat.objects)] do
         if k in cat.objects[x] then
              return(x);
         fi;
    od;
end );

InstallGlobalFunction(ConcreteCategory, function(arg)
    local output, domains, codomains, x, m, cod, dom, included, obj, nums;
    output:=rec(generators:=arg[1], operations:=ConcreteCategoryOps);
    if Length(arg)=2 then
        #Checks if object sets overlap or entries repeated (catches ([1,2],[2,3,4]) or ([1,2,2]))
        nums:=[];
        for x in arg[2] do
            for m in x do
                if not (m in nums) then
                    Add(nums,m);
                else
                    Error("One or more entries is duplicated within one or more objects.");
                fi;
            od;
        od;

        #Computes domains/codomains of morphisms
        domains:=[];
        for m in arg[1] do
            Add(domains,SupportOfMorphism(m));
        od;
        domains:=SSortedList(domains);
        codomains:=[];
        for x in domains do
            for m in output.generators do
                if SupportOfMorphism(m)=x then
                    Add(codomains, List(x,a->m[a]));
                fi;
            od;
        od;
        codomains:=SSortedList(codomains);
        
        #Ensures domains/codomains of morphisms in provided objects list
        for dom in domains do
            included:=false;
            for obj in arg[2] do
                if IsSubset(obj,dom) then
                    included:=true;
                    break;
                fi;
            od;
            if included=false then
                Error("One or more morphisms have domains not included in objects list.");
            fi;
        od;
        for cod in codomains do
            included:=false;
            for obj in arg[2] do
                if IsSubset(obj,cod) then
                    included:=true;
                    break;
                fi;
            od;
            if included=false then
                Error("One or more morphisms have codomains not included in objects list.");
            fi;
        od;

        output.objects:=SSortedList(arg[2]);
    
    elif Length(arg)=1 then
        Objects(output);
    fi;

    Domains(output);
    return(output);
end );

##############################################################################
##
## Morphisms(cat) returns an l x l matrix, where l is the number of objects
## in the category cat, and where the i,j entry is a list of the
## morphisms from object i to
## object j.
##
##############################################################################
InstallGlobalFunction( Morphisms, function(cat)
    local n, genmat, g, mormat, oldlength, newlength, i, j, k, h, templist;
    if IsBound(cat.morphisms) then return(cat.morphisms);
    fi;
    if not IsBound(cat.objects) then Objects(cat);
    fi;
    n:=Length(cat.objects);
    genmat:=EmptyMat(n,n);
    mormat:=EmptyMat(n,n);
    for g in cat.generators do
        Add(genmat[Origin(cat,g)][Terminus(cat,g)],g);
    od;
    for i in [1..n] do
        Add(mormat[i][i],IdentityMorphism(cat.objects[i]));
    od;
    oldlength:=0;
    newlength:=Sum(List(mormat,x->Sum(List(x,Length))));
    while oldlength < newlength do
        oldlength:=newlength;
        for i in [1..n] do
            for j in [1..n] do
                templist:=[];
                for k in [1..n] do
                    for g in genmat[k][j] do
                        for h in mormat[i][k] do
                            Add(templist, Composition(h,g));
                        od;
                    od;
                 od;
                 Append(mormat[i][j],templist);
                 mormat[i][j]:=SSortedList(mormat[i][j]);
            od;
        od;
        newlength:=Sum(List(mormat,x->Sum(List(x,Length))));
    od;
    cat.morphisms:=mormat;
    return(mormat);
end);

InstallGlobalFunction( EndomorphismGroups, function( cat )
    cat.endomorphismgroups := List( cat.objects, x -> EndomorphismGroup( cat, x ) );
end );

InstallGlobalFunction( FI, function(n)
    local objectlist, morphismlist, i, j, x, m;
    morphismlist:=[];
    
    objectlist:=[];
    j:=0;
    for i in [1..n+1] do
        Add(objectlist, [j+1..j+i]);
        j:=j+i;
    od;
    
    for x in objectlist do
        if Length(x)<=n then
            m:=[];
            for i in x do
                m[i]:=i+Length(x);
            od;
            Add(morphismlist,m);
        fi;
        
        if Length(x)>=3 then
            m:=[];
            for i in x do
                m[i]:=i;
            od;
            m[x[2]]:=x[2]+1;
            m[x[3]]:=x[3]-1;
            Add(morphismlist,m);
        fi;
        
        if Length(x)>=4 then
            m:=[];
            m[x[1]]:=x[1];
            m[x[Length(x)]]:=x[2];
            for i in [x[2]..x[Length(x)-1]] do
                m[i]:=i+1;
            od;
            Add(morphismlist,m);
        fi;
    od;
    
    return ConcreteCategory(morphismlist, objectlist);
end );

InstallGlobalFunction( FI2, function(n)
    local objectlist, morphismlist, i, j, k, m;
    morphismlist:=[];
    objectlist:=[];
    j:=0;
    for i in [1..n+1] do
        Add(objectlist, [j+1..j+i]);
        
        if i<=n then
            m:=[];
            for k in [j+1..j+i] do
                m[k]:=k+i;
            od;
            Add(morphismlist,m);
        fi;
        
        if i>=3 then
            m:=[];
            for k in [j+1..j+i] do
                m[k]:=k;
            od;
            m[j+2]:=j+3;
            m[j+3]:=j+2;
            Add(morphismlist,m);
        fi;
        
        if i>=4 then
            m:=[];
            m[j+1]:=j+1;
            m[j+i]:=j+2;
            for k in [j+2..j+i-1] do
                m[k]:=k+1;
            od;
            Add(morphismlist,m);
        fi;
        j:=j+i;
    od;
    
    return ConcreteCategory(morphismlist, objectlist);
end );


