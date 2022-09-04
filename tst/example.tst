gap> # We first set up a category with 2 objects and 9 morphisms:
gap> c3c3:=ConcreteCategory([[2,3,1],[4,5,6],[,,,5,6,4]]);;
gap> # The following constructs a representation:
gap> one:=One(GF(3));;
gap> d:=[[1,1,0,0,0],[0,1,1,0,0],[0,0,1,0,0],[0,0,0,1,1],[0,0,0,0,1]]*one;;
gap> e:=[[0,1,0,0],[0,0,1,0],[0,0,0,0],[0,1,0,1],[0,0,1,0]]*one;;
gap> f:=[[1,1,0,0],[0,1,1,0],[0,0,1,0],[0,0,0,1]]*one;;
gap> nine:=CatRep(c3c3,[d,e,f],GF(3));;
gap> # The dimension vector shows that the space
gap> # associated to object 1 is 5-dimensional,
gap> # and the space associated to object 2 is 4-dimensional.
gap> # The representation is indecomposable:
gap> Length(Decompose(nine));
1
gap> fortyone:=TensorProductRep(nine,nine);;
gap> # The output from the last command was suppressed,
gap> # but the representation it produces has a space of
gap> # dimension 25 at object 1, and a space of dimension
gap> # 16 at object 2. We start to examine it by finding
gap> # its indecomposable summands.
gap> d:=Decompose(fortyone);;
gap> List(d,x->List(x,Length));
[ [ 3, 0 ], [ 3, 1 ], [ 3, 3 ], [ 3, 3 ], [ 0, 3 ], [ 3, 0 ], [ 3, 0 ], 
  [ 3, 0 ], [ 1, 3 ], [ 3, 3 ] ]
gap> # Thus fortyone has 10 indecomposable summands,
gap> # with dimension lists as above. We examine the}
gap> # structure of the third summand.
gap> six:=SubmoduleRep(fortyone,d[3]);;
gap> SumOfImages(ConstantRep(c3c3,GF(3)),six);
[ [ [ 0*Z(3), 0*Z(3), Z(3)^0 ] ], [  ] ]
gap> # The representation six consists of a module homomorphism
gap> # for a group of order 3 acting on two spaces of dimension 3.
gap> # The last calculation shows that it is not an isomorphism
gap> # between two modules. We do the same test to the Yoneda
gap> # representation (= representable functor) associated to object 1.
gap> y:=YonedaRep(c3c3,1,GF(3));;
gap> SumOfImages(ConstantRep(c3c3,GF(3)),y);
[ [ [ Z(3)^0, Z(3)^0, Z(3)^0 ] ], [ [ Z(3)^0, Z(3)^0, Z(3)^0 ] ] ]
gap> # This shows that six is not isomorphic to the Yoneda
gap> # representation. We take the opportunity to construct
gap> # the simple representation associated to object 1.
gap> s1:=SubmoduleRep(six,[ [ [ 0*Z(3), 0*Z(3), Z(3)^0 ] ], [ ] ]);;
gap> five:=QuotientRep(six,[ [ [ 0*Z(3), 0*Z(3), Z(3)^0 ] ], [ ] ]);;
gap> # The next calculation shows that the 3-dimensional representation
gap> # of C_3 associated to object 1 is a single copy of the regular
gap> # representation of C_3.
gap> SumOfImages(s1,six);
[ [ [ 0*Z(3), 0*Z(3), Z(3)^0 ] ], [  ] ]
gap> # The next calculation shows that the quotient representation
gap> # five maps its module at object 1 monomorphically to the module
gap> # at object 2, which must either be indecomposable of dimension
gap> # 3, or else the direct sum of indecomposables of dimension 2 and 1.
gap> SumOfImages(s1,five);
[ [  ], [  ] ]
gap> SumOfImages(ConstantRep(c3c3,GF(3)),five);
[ [ [ 0*Z(3), Z(3)^0 ] ], [ [ 0*Z(3), 0*Z(3), Z(3)^0 ] ] ]
gap> # The final calculation shows that the module at object
gap> # 2 for six is indecomposable of dimension 3. We now have
gap> # sufficient information to describe six completely.
gap> SumOfImages(six,ConstantRep(c3c3,GF(3)));
[ [  ], [ [ Z(3)^0 ] ] ]
