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
#! @Chapter Introduction

#! @Section Overview of the package reps
#! 
#! The package `reps` described here sets up a framework
#! within which matrix representations of finite groups,
#! primarily in positive characteristic, can be stored and
#! examined. Sufficient information is stored so that calculations
#! can be performed, and at the same time little more than this
#! information is stored. It is important that the user should
#! be able to understand quickly the syntax of commands and what
#! the software is doing, and also that the software should not
#! depend unduly on constructions in other parts of &GAP; that may
#! change.
#! 
#! The basic object of computation is a record with a name such
#! as `rep`, and fields which include `rep.group`, `rep.genimages`,
#! `rep.field` and `rep.dimension`. The first of these is the group
#! being represented, and the second is a list of matrices over
#! `rep.field` which represent the action of the group elements in
#! the list `GeneratorsOfGroup(rep.group)`. Since many of the
#! algorithms used and some of the properties investigated
#! (such as projectivity) depend on the group being represented
#! it is important to store the group. By contrast, the implementation
#! of the meataxe already present in &GAP; does not store the group.
#! With these record fields it is also possible to distinguish
#! between the zero representation of an arbitrary group and the
#! different representations of the identity group. In the present
#! implementation certain commands will only work when `rep.group` is
#! a permutation group.
#! 
#! The intention is to make available all algorithms for
#! handling group representations within the above framework.
#! These include the basic meataxe operations for finding a proper
#! invariant subspace and Norton's irreducibility test, but also
#! include many other algorithms, many of them based on the operation
#! of taking fixed points. It is important to understand the
#! limitations of the various algorithms in using them effectively.
#! Thus the meataxe is at its best when given a representation with
#! few, non-isomorphic, composition factors. To find the socle of a
#! representation of a $p$-group in characteristic $p$, it may be better
#! to use `FixedPoints(rep);` which solves some linear equations,
#! rather than `SocleRep(rep);` which calls the meataxe. The routines
#! which work with spaces of homomorphisms, such as Decompose, tend
#! to be limited by the dimension of the representations; and so on.
#!
#! The fastest way to become familiar with this package is to
#! read through the worked examples of its use which appear after
#! the list of available commands.
#!
#!
#! @Subsection Installation
#!
#! To install this package, refer to the installation instructions in
#! the reps file in the source code. It is located here:
#! <URL>https://www-users.cse.umn.edu/~webb/GAPfiles/reps</URL>.


#! @Chapter Commands that return a representation

#! @Chapter Commands that return lists of vectors

#! @Chapter Commands that handle output

#! @Chapter Commands for homomorphisms

#! @Chapter Trace, norm and principal idempotent

#! @Chapter Test

#! @Chapter Interface to the Meataxe