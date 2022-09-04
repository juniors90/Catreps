# ***************************************************************************
#                                   Reps
#       Copyright (C) 2005 Peter Webb 
#       Copyright (C) 2006 Peter Webb, Robert Hank, Bryan Simpkins 
#       Copyright (C) 2007 Peter Webb, Dan Christensen, Brad Froehle
#       Copyright (C) 2020 Peter Webb, Moriah Elkin
#       Copyright (C) 2022 Peter Webb, Ferreira Juan David
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#
#  The overall structure of the reps package was designed and most of it
#  written by Peter Webb <webb@math.umn.edu>, who is also the maintainer. 
#  Contributions were made by
#  Dan Christensen, Roland Loetscher, Robert Hank, Bryan Simpkins,
#  Brad Froehle and others.
# ***************************************************************************
#
# LoadPackage("Catreps", "0", false);
#! @BeginChunk Example_1
#! @BeginExample
LoadPackage("Catreps", "0", false);
true
# We first set up a category with 2 objects and 9 morphisms:
c3c3:=ConcreteCategory([[2,3,1],[4,5,6],[,,,5,6,4]]);;
# The following constructs a representation:
one:=One(GF(3));;
d:=[[1,1,0,0,0],[0,1,1,0,0],[0,0,1,0,0],[0,0,0,1,1],[0,0,0,0,1]]*one;;
e:=[[0,1,0,0],[0,0,1,0],[0,0,0,0],[0,1,0,1],[0,0,1,0]]*one;;
f:=[[1,1,0,0],[0,1,1,0],[0,0,1,0],[0,0,0,1]]*one;;
nine:=CatRep(c3c3,[d,e,f],GF(3));;
# The dimension vector shows that the space
# associated to object 1 is 5-dimensional,
# and the space associated to object 2 is 4-dimensional.
# The representation is indecomposable:
Length(Decompose(nine));
#! 1
fortyone:=TensorProductRep(nine,nine);;
# The output from the last command was suppressed,
# but the representation it produces has a space of
# dimension 25 at object 1, and a space of dimension
# 16 at object 2. We start to examine it by finding
# its indecomposable summands.
d:=Decompose(fortyone);;
List(d,x->List(x,Length));
#! [ [ 3, 0 ], [ 3, 1 ], [ 3, 3 ], [ 3, 3 ], [ 0, 3 ],
#!   [ 3, 0 ], [ 3, 0 ], [ 3, 0 ], [ 1, 3 ], [ 3, 3 ] ]
# Thus fortyone has 10 indecomposable summands,
# with dimension lists as above. We examine the}
# structure of the third summand.
six:=SubmoduleRep(fortyone,d[3]);;
SumOfImages(ConstantRep(c3c3,GF(3)),six);
#! [ [ [ 0*Z(3), 0*Z(3), Z(3)^0 ] ], [  ] ]
# The representation six consists of a module homomorphism
# for a group of order 3 acting on two spaces of dimension 3.
# The last calculation shows that it is not an isomorphism
# between two modules. We do the same test to the Yoneda
# representation (= representable functor) associated to object 1.
y:=YonedaRep(c3c3,1,GF(3));;
SumOfImages(ConstantRep(c3c3,GF(3)),y);
#! [ [ [ Z(3)^0, Z(3)^0, Z(3)^0 ] ], [ [ Z(3)^0, Z(3)^0, Z(3)^0 ] ] ]
# This shows that six is not isomorphic to the Yoneda
# representation. We take the opportunity to construct
# the simple representation associated to object 1.
s1:=SubmoduleRep(six,[ [ [ 0*Z(3), 0*Z(3), Z(3)^0 ] ], [ ] ]);;
five:=QuotientRep(six,[ [ [ 0*Z(3), 0*Z(3), Z(3)^0 ] ], [ ] ]);;
#! The next calculation shows that the 3-dimensional representation
#! of C_3 associated to object 1 is a single copy of the regular
#! representation of C_3.
SumOfImages(s1,six);
#! [ [ [ 0*Z(3), 0*Z(3), Z(3)^0 ] ], [  ] ]
# The next calculation shows that the quotient representation
# five maps its module at object 1 monomorphically to the module
# at object 2, which must either be indecomposable of dimension
# 3, or else the direct sum of indecomposables of dimension 2 and 1.
SumOfImages(s1,five);
#! [ [  ], [  ] ]
SumOfImages(ConstantRep(c3c3,GF(3)),five);
#! [ [ [ 0*Z(3), Z(3)^0 ] ], [ [ 0*Z(3), 0*Z(3), Z(3)^0 ] ] ]
# The final calculation shows that the module at object
# 2 for six is indecomposable of dimension 3. We now have
# sufficient information to describe six completely.
SumOfImages(six,ConstantRep(c3c3,GF(3)));
#! [ [  ], [ [ Z(3)^0 ] ] ]
#! @EndExample
#! @EndChunk