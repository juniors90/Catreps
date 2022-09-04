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

InstallGlobalFunction( DimHom, function(M, N )
    return Length( HomBasis( M, N ) );
end );

InstallGlobalFunction( HomBasis, function(M,N)
    if not IsBound(M.operations.HomBasis) and not IsBound(N.operations.HomBasis) then
        Error("for usage, see ?HomBasis");
	fi;
    return M.operations.HomBasis( M, N );
end );

InstallGlobalFunction( CatHomBasis, function( M, N )
    local dimM, dimN, dM, dN, dMN, dimsumM, dimsumN, dimsumMN, r, s, i, j, k, l, v, f, domain, codomain;
    # dimM := M.dimension;
    # dimN := N.dimension;
    dimsumM:=[];  dimsumM[1]:=0;
    dimsumN:=[];  dimsumN[1]:=0;
    dimsumMN:=[]; dimsumMN[1]:=0;
    domain:=M.category.domain; codomain:=M.category.codomain;
    r := Length(M.genimages);
    dM:=0; dN:=0; dMN:=0;
    for i in [1..Length(M.dimension)] do
        dM:=dM+M.dimension[i];
        dN:=dN+N.dimension[i];
        dMN:=dMN+M.dimension[i]*N.dimension[i];
        Add(dimsumM,dM);
        Add(dimsumN,dN);
        Add(dimsumMN,dMN);
    od;
    dimM := dimsumM[Length(M.dimension)+1];
    dimN := dimsumN[Length(M.dimension)+1];
    v:=NullMat(dimsumMN[Length(M.dimension)+1], dimM*dimN*r, M.field);
    for l in [1..r] do
        for s in [1..N.dimension[codomain[l]]] do
            for j in [1..N.dimension[domain[l]]] do
                for i in [1..M.dimension[domain[l]]] do
                    v[dimsumMN[domain[l]]+(i-1)*N.dimension[domain[l]]+j]
                    [(l-1)*dimM*dimN+dimsumM[domain[l]]*dimN+M.dimension[domain[l]]*dimsumN[codomain[l]]+(i-1)*N.dimension[codomain[l]]+s]
                    :=v[dimsumMN[domain[l]]+(i-1)*N.dimension[domain[l]]+j]
                    [(l-1)*dimM*dimN+dimsumM[domain[l]]*dimN+M.dimension[domain[l]]*dimsumN[codomain[l]]+(i-1)*N.dimension[codomain[l]]+s]
                    - N.genimages[l][j][s];
                od;
            od;
        od;
        for r in [1..M.dimension[domain[l]]] do
            for j in [1..N.dimension[codomain[l]]] do
                for i in [1..M.dimension[codomain[l]]] do
                    v[dimsumMN[codomain[l]]+(i-1)*N.dimension[codomain[l]]+j]
                    [(l-1)*dimM*dimN+dimsumM[domain[l]]*dimN+M.dimension[domain[l]]*dimsumN[codomain[l]]+(r-1)*N.dimension[codomain[l]]+j]
                    :=v[dimsumMN[codomain[l]]+(i-1)*N.dimension[codomain[l]]+j]
                    [(l-1)*dimM*dimN+dimsumM[domain[l]]*dimN+M.dimension[domain[l]]*dimsumN[codomain[l]]+(r-1)*N.dimension[codomain[l]]+j]
                    + M.genimages[l][r][i];
                od;
            od;
        od;
    od;
    f := SafeNullspaceMat( v, M.field );
    return( List( f, x -> ExtractHom( x, M.dimension, N.dimension ) ) );
end );
