# Commands for categories:

#! Rep( group, matrices, optional field)
#! This function creates a representation provided the group is not the
#! identity group and the representation is not the zero representation.
#! Modification by Craig Corsi and Peter Webb (University of Minnesota)
#! 2016 to allow the optional field argument.
DeclareGlobalFunction( "Rep" );

#! SupportOfMorphism(m) returns the list of positions at which the 
#! list m is defined.
#!
#! SupportOfMorphism(mapping) returns the list of positions
#! at which the list 'mapping' is defined.
DeclareGlobalFunction( "SupportOfMorphism" );

#! IdentityMorphism(l) returns the identity morphism on the set l.
#! A morphism is stored as a list of the images of its values on a set
#! of numbers, which form its domain of definition, and are taken to be
#! an object in a concrete category. At elements of other objects the
#! morphism will be undefined.
#!
#! IdentityMorphism(l) returns the identity morphism on the set l.
DeclareGlobalFunction( "IdentityMorphism" );

#! Composition(f,g) returns the composition of the
#! functions f and g, expressed as lists of their values.
#!
#! Composition(f,g) returns the composition of the
#! functions f and g, expressed as lists of their values.
DeclareGlobalFunction( "Composition" );

#! IsComposable(f, g) returns true if the functions
#! f and g, expressed as lists of their values, can
#! be composed and false otherwise.
#!
#! IsComposable(f, g) returns true if the functions
#! f and g, expressed as lists of their values, can
#! be composed and false otherwise.

DeclareGlobalFunction( "IsComposable" );

#! Objects(cat) returns the objects of the concrete
#! category cat, as a list of sets.
#!
#! Objects(cat) returns the objects of the concrete category
#! cat, as a list of sets. At the moment it will not work
#! unless for every object there is at least one generator
#! morphism whose support is that object.
DeclareGlobalFunction( "Objects" );

#! Origin(cat,m) returns the position in cat.objects
#! of the domain of the morphism m.
#!
#! Origin( cat, m ) returns the position in cat.objects of the
#! domain of the morphism m.
DeclareGlobalFunction( "Origin" );

#! Terminus( cat, m ) returns the position in cat.objects
#! of the domain of the morphism m.
#!
#! Terminus( cat, m ) returns the position in cat.objects of
#! the domain of the morphism m.
DeclareGlobalFunction( "Terminus" );

#! ConcreteCategory(list of mappings (, list of sets)). . sets up
#! the record of the category which has the list of mappings as
#! its generator morphisms, with the optional second argument
#! as its list of objects. If the second argument is missing,
#! the objects are taken to be the domains of the morphisms.
#!
#! ConcreteCategory(list of functions, (list of sets))
#! There are optionally one or two arguments. The first is a
#! list of generating functions of the category, the second
#! is a list of the objects of the category.
#!
#! The function starts a record for a concrete category.
#! If there is only one argument, the objects are taken
#! to be the domains of the generator morphisms, so for
#! every object there should be at least one generator morphism
#! whose support is that object.
#!
#! It could be the identity morphisms, but doesn't have to be.
#!
#! Written by Peter Webb 2008, Moriah Elkin 2018.
DeclareGlobalFunction( "ConcreteCategory" );

#! Morphisms( cat ) returns an l x l matrix, where l is
#! the number of objects in the category cat, and where
#! the i,j entry is a list of the morphisms from object
#! i to object j.
#!
#! Morphisms( cat ) returns an l x l matrix, where l is
#! the number of objects in the category cat, and where
#! the i,j entry is a list of the morphisms from object
#! i to object j.
DeclareGlobalFunction( "Morphisms" );

#! EndomorphismGroups(cat) creates a field of the record
#! cat listing the endomorphism groups of objects. It is
#! assumed these are groups, and that their generators
#! appear among the generators of the category.
DeclareGlobalFunction( "EndomorphismGroups" );

#! FI(n). .returns the full subcategory of the category
#! FI of finite sets and injective maps, whose objects
#! are 0,...,n. The empty set 0 is represented as a
#! list [1], and the unique morphism from the empty set
#! to each object is represented as a map sending 1 to
#! the first element in that object, which is otherwise
#! ignored. Thus the first element of each object is a
#! place holder for the empty set.

#!
#! FI(n) and FI2(n) interchangeably return a record
#! for the category FI with objects 0...n. O is represented
#! by the first object ( [1] ), and its morphisms
#! correspond to the first element in every object,
#! which is otherwise ignored (first elements only map to
#! first elements). The category FI is the category of
#! finite sets with injective maps that has featured in
#! the theory of representation stability.
#!
#! Written by Moriah Elkin July 2018.
DeclareGlobalFunction( "FI" );

# ------------------------------------------------------------
DeclareGlobalFunction( "FI2" );
