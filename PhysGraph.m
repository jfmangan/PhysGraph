(* ::Package:: *)

(* ::Subsection:: *)
(*Begin package*)


BeginPackage["PhysGraph`",{"Amplitudes`"}]


(* ::Subsection:: *)
(*Load packages and general notes*)


(* ::Text:: *)
(*VerGraph means a vertex graph.*)
(**)
(*VerF[a,b,c...] is the Head used in VerGraph for representing an oriented n-pt vertex.  Cyclical order doesn't matter but swaps are associated with a minus sign just like f^abc.  So the name is "VerF" where "Ver" is for vertex and "F" is for fabc.  Note that external/single vertices like VerF[p[1]] are included in VerGraph.*)
(*VerGraph[{VerF[...],VerF[...]...}] represents a graph with vertices VerF[].  It is not an MMA graph.*)
(**)
(*^^This is all true for VerF[a,b,c] but for 4 or more legs some of this isn't true.  This can become an issue for cuts where internal blobs have 4 or more legs but this issue is handled in GenInternalLabelingsOfGraphOfCuts.*)
(**)
(*I'm choosing a clockwise convention for the edges a,b,c... in VerF[a,b,c...].  MMA internally choose counterclockwise I think.*)
(**)
(*All outgoing convention.*)
(**)
(*Label edges in vertices with signs for outgoing (positive) or incoming (negative).  This means that external vertices look like VerF[-a[x]] if you follow the all outgoing convention.  This will get messed up in a cut where you'll have a pair like VerF[-a[x]] and VerF[a[x]].*)
(**)
(*If you have graphs like VerGraph[{VerF[-a[x]...}] or VerGraph[{VerF[-p[x]...}] then you need to add a or p to LVecHeads!*)


(* ::Text:: *)
(*Remember that graphs store an internal sign based on the orientation of the legs in all of the VerFs.  So when you display a graph there could still be a hidden sign compared to what the graph looks like.*)


(* ::Text:: *)
(*Things to always worry about:*)
(**signs, especially on graph isomorphisms*)
(**symmetry factors*)


(* ::Text:: *)
(*Three kinds of functional relations:*)
(*1) From graph syms/automorphisms num[1][p1,p2...] == num[1][p2,p7...]*)
(*2) Jacobi relations num1+num1+num2=0*)
(*3) Unitarity cuts num[1][...]+num[2][...]=A4()A4()A3()*)


(* ::Subsection:: *)
(*Usage notes*)


VerGraph::usage="VerGraph[VerList] is the Head associated with a vertex graph of vertices {VerF[a,b,c],VerF[-c,d,e]...}";

VerF::usage="VerF[a,b,c] is the Head associated with a vertex with momentum/edge labels {a,b,c}.  The 'F' in 'VerF' stands for f^abc because (at least at 3pt) the arguments of VerF are treated as though they are cyclically invariant and totally antisymmetric.  For quartic (or higher) vertex, the labels are still treated as cyclically invariant but no other symmetries are assumed.  This is important when working with higher point blobs in cuts.  This is handled by GenInternalLabelingsOfGraphOfCuts.";

VerGraphToMMAGraph::usage="VerGraphToMMAGraph[VerGraph] converts a VerGraph to MMA's default graph type.";

MultigraphDisambiguate::usage="MultigraphDisambiguate[MMAGraph] takes a MMA Graph and adds dots to edges (spurious vertices to propagators) to turn a Mulitgraph into a normal graph.";

VerGraphHashCode::usage="VerGraphHashCode[VerGraph] returns a canonicalized MMA graph.  This may involved dotting edges.  This function is memoized but the cache is cleared by RunGraphSetup.";

GetExtVerts::usage="GetExtVerts[VerGraph] returns the external vertices (or equivalently the outgoing momenta) of a VerGraph.";

GenCubicTopologies::usage="GenCubicTopologies[n,l,MomHead] will generate the n-pt l-loop distinct cubic graph topologies with momenta labeled by MomHead.";

MomConsRepRules::usage="MomConsRepRules[VerGraph,MomToSolveFor:{}] gives momentum conservation rules for VerGraph where the (optional) list of momenta MomToSolveFor have been eliminated.";

dRepRules::usage="dRepRules[VerGraph,AdditionalConstraints:{},MomToSolveFor:{}] returns the replacement rules on the dot products of momenta d[p[x],p[y]] that enforces on-shell constraints.  The replacement rule does not work on polarization vectors.  Note that you need to call MomConsRepRules first.  AdditionalConstraints is an optional list of constraints to impose.  This is useful for doing cuts.  MomToSolveFor is an optional list of momenta to explicitly eliminate.  Since you need to call MomConsRepRules first, make sure to use the same MomToSolveFor in both functions or you will get inconsistent bases from MomConsRepRules and dRepRules.";

RepGraphKin::usage="RepGraphKin[graph][expr] conserves momentum and replaces the dot products of momenta associated with graph in the expression expr.  This is equivalent to a call to MomConsRepRules and then dRepRules.";

MandelstamBasisOfGraph::usage="MandelstamBasisOfGraph[VerGraph] returns the basis of momentum dot products associated with the graph VerGraph.  This is the basis used by MomConsRepRules, dRepRules, and RepGraphKin.";

TadpoleQ::usage="TadpoleQ[VerGraph] returns True if there is a tadpole in the graph VerGraph";

TadpoleOrBubbleOnExtLegQ::usage="TadpoleOrBubbleOnExtLegQ[VerGraph] returns True if the graph VerGraph contains either a tadpole or a bubble on an external leg.";

GraphSignature::usage="GraphSignature[VerGraph] returns the signature (+1 or -1) of the graph VerGraph.  This sign has to do with how vertices of VerGraph are oriented.";

VerGraphToMMAGraph::usage="VerGraphToMMAGraph[VerGraph] converts a VerGraph to an MMA graph.";

VerGraphToEdgeGraph::usage="VerGraphToEdgeGraph[VerGraph] converts a VerGraph to an EdgeGraph.";

EdgeGraphToVerGraph::usage="EdgeGraphToVerGraph[EdgeGraph] converts an EdgeGraph to a VerGraph";

PrintEdgeGraph::usage="PrintEdgeGraph[EdgeGraph] prints the edge graph EdgeGraph along with momentum/edge labels and their directions.  The output is far from perfect.";

VerGraphMultigraphDisambiguate::usage="VerGraphMultigraphDisambiguate[VerGraph] adds dots (fictitious 2pt vertices) to edges so that VerGraph isn't a Multigraph.  The function returns a VerGraph.";

FindVerGraphIsomorphisms::usage="FindVerGraphIsomorphisms[VerGraphSource,VerGraphTarget] finds an isomorphism from VerGraphSource to VerGraphTarget.";

FindVerGraphSpanningAutomorphisms::usage="FindVerGraphSpanningAutomorphisms[VerGraph] finds a set of automorphisms on VerGraph.  I try to find the smallest set of automorphisms so that when you compose them you can get any automorphism however there is still some redundancy.  The reason for all of this is to try to generate as few redundant constraints when imposing graph symmetries on numerator ansatze.";

RunGraphSetup::usage="RunGraphSetup[n,l,MomentumHead] applies AddLVecHeads to MomentumHead and populates CubicBasisGraphs, CanonicalNumeratorVariables, BasisGraphNumeratorPairs, and BasisMMAGraphNumeratorPairs.  This function must be run before generating graph symmetries, Jacobi relations, or doing cuts.  CubicBasisGraphs is the set of cubic n-pt l-loop graphs with no tadpoles or bubbles on external legs.  CanonicalNumeratorVariables is a list that loos like {p[1], p[2]...} where p is MomentumHead.  The numerator associated with CubicBasisGraphs[[i]] is num[i][CanonicalNumeratorVariables].  BasisGraphNumeratorPairs is the List of {CubicBasisGraph[[i]],num[i][CanonicalNumeratorVariables]}.  BasisMMAGraphNumeratorPairs is the same as BasisGraphNumeratorPairs but each CubicBasisGraph replaced with its VerGraphHashCode.  RunGraphSetup also clears the caches of the memoized functions GenNkMCGraphsOfCuts and VerGraphHashCode.  If you change CubicBasisGraphs then you definitely need to clear the cache for GenNkMCGraphsOfCuts and you probably don't need the cached values of VerGraphHashCode any more either.";

CubicBasisGraphs::usage="CubicBasisGraphs is the set of cubic n-pt l-loop graphs with no tadpoles or bubbles on external legs.  CubicBasisGraphs is populated by RunGraphSetup.";

CanonicalNumeratorVariables::usage="CanonicalNumeratorVariables is a list that loos like {p[1], p[2]...} where p is MomentumHead.  The numerator associated with CubicBasisGraphs[[i]] is num[i][CanonicalNumeratorVariables].  CanonicalNumeratorVariables is populated by RunGraphSetup.";

BasisGraphNumeratorPairs::usage="BasisGraphNumeratorPairs is the List of {CubicBasisGraph[[i]],num[i][CanonicalNumeratorVariables]}.  BasisGraphNumeratorPairs is populated by RunGraphSetup.";

BasisMMAGraphNumeratorPairs::usage="BasisMMAGraphNumeratorPairs is the same as BasisGraphNumeratorPairs but each CubicBasisGraph replaced with its VerGraphHashCode and each num[i][...] replaced with just num[i].  BasisMMAGraphNumeratorPairs is populated by RunGraphSetup.";

num::usage="num[4][{p[2], p[5]...}] is the Head used to represent a graph numerator.  The argument is a list of momentum labels.  num[i][CanonicalNumeratorLabels] corresponds to the numerator of CubicBasisGraph[[i]].  Relabeling the momentum labels of num[i] is supposed to correspond to an isomorphism (or momentum relabeling) of CubicBasisGraph[[i]].";

PinchSpecifiedEdge::usage="PinchSpecifiedEdge[VerGraph,edge] pinches edge (for example edge=p[7]) in VerGraph.  In other words collapses this function collapses propagators.";

PinchOneEdge::usage="PinchOneEdge[VerGraph] returns every graph that can be obtained by pinching one edge (or collpasing one propagator) of VerGraph.";

PinchOneEdgeOnSetOfGraphs::usage="PinchOneEdgeOnSetOfGraphs[VerGraphList] takes a list of VerGraphs and returns every graph that can be obtained by pinching one edge (or collapsing one propagator) of one of the input graphs.  This is equivalent to running PinchOneEdge on all of the graphs and modding out by isomorphisms.";

BlowUpVertexIntoJacobiTriplet::usage="BlowUpVertexIntoJacobiTriplet[VerGraph] takes a VerGraph that is assumed to have exactly one quartic vertex.  It blows up that vertex in the s, t, and u channels and returns the resulting three VerGraphs.  The leg inserted when blowing up the graph is EdgeInsertedInJacobi which is set by RunJacobiSetup.";

FindVerGraphNumerator::usage="FindVerGraphNumerator[VerGraph] takes any labeling of VerGraph and returns the associated numerator num[i][labeling] with the corresponding labeling.";

FindConstructibleNumerators::usage="FindConstructibleNumerators[nums] takes a list of numerators in simplified notation {num[1],num[4]...} and generates every numerator that can be obtained by using Jacobi relations.  This function relies on SimplifiedJacobis which is set up by RunJacobiSetup.";


RunJacobiSetup::usage="RunJacobiSetup[IncludeTadpoleBELJacobis] runs all of the Jacobi setup code including setting up EdgeInsertedInJacobi, JacobiRelationsOnGraphs, JacobiRelationsOnNumerators, and SimplifiedJacobis. If IncludeTadpoleBELJacobis is False then you only keep Jacobis of length 3, that is, you ignore solving any Jacobi involving a tadpole or a bubble on external leg (BEL). If IncludeTadpoleBELJacobis is True then you solve all Jacobi relations (even those involving tadpoles or BELs) but you hardcode the numerators of tadpoles and BELs to zero.";


EdgeInsertedInJacobi::usage="Is the edge (like p[17]) that is inserted into a graph when blowing up a quartic vertex to obtain the Jacobi relations on a graph.  EdgeInsertedInJacobi is used in BlowUpVertexIntoJacobiTriplet and is set up by RunJacobiSetup.";

JacobiRelationsOnGraphs::usage="JacobiRelationsOnGraphs represents the Jacobi relations at the VerGraph level.  JacobiRelationsOnGraphs is a table of Jacobi relations where each element is a list of the graphs contributing to that Jacobi relations.  JacobiRelationsOnGraphs is populated by RunJacobiSetup where you get to decide if you want to keep/solve Jacobis involving tadpole and bubble on external leg (BEL) graphs.";

JacobiRelationsOnNumerators::usage="JacobiRelationsOnNumerators is the set of functional Jacobi relations at the numerator (num[...][...]) level.  Momentum conservation has been used to eliminate EdgeInsertedInJacobi.  JacobiRelationsOnNumerators is a list of the numerators contributing to each Jacobi relation.  If the Jacobi relation is satisfied then the Total of that list should be zero.  JacobiRelationsOnNumerators is set up by RunJacobiSetup.";

SimplifiedJacobis::usage="SimplifiedJacobis is used in finding a basis of numerators.  SimplifiedJacobis is the list of Jacobi relations in a simplified notation {num[i],num[j],num[j]} or {num[i],num[j]}.  Self Jacobis like {num[i],num[i],num[i]} or {num[i],num[i]} are not kept because they are not helpful in finding a basis of numerators.  SimplifiedJacobis is set up by RunJacobiSetup.";


FindConstructibleNumerators::usage="FindConstructibleNumerators[nums] takes a list of numerators in simplified form {num[1], num[7]...} and returns every numerator that can be obtained from nums by applying Jacobi relations.  FindConstructibleNumerators relies on SimplifiedJacobis which is set up by RunJacobiSetup.";

FindSomeBasisOfNumerators::usage="FindSomeBasisOfNumerators populates nBasis by finding a set of basis/master numerators. You must RunJacobiSetup before using FindSomeBasisOfNumerators.";

nBasis::usage="nBasis is a set of numerators in simplified notation {num[1], num[7]...}.  Every numerator can be obtained through Jacobi relations on nBasis.  You can populate nBasis by hand or use FindSomeBasisOfNumerators to generate the basis.  Note that nBasis is not unique in general.  You should RunJacobiSetup before setting up nBasis.  If you construct nBasis by hand you can check that it produces a valid basis by making sure that FindConstructibleNumerators[nBasis] gives all numerators.  You must RunJacobiSetup before using FindConstructibleNumerators.";

SimplifiedNumsToVerGraphs::usage="SimplifiedNumsToVerGraphs[NumList] takes a list of numerators in simplified notation {num[1], num[7]...} and returns the associated VerGraphs.  For example SimplifiedNumsToVerGraphs[nBasis] is helpful in visually checking that you have the right basis graphs.";

NumeratorsToMasterNumerators::usage="NumeratorsToMasterNumerators[expr] replaces every num[...][...] in expr with the appropriately relabeled num from nBasis.  You must RunNumeratorsToMasterNumeratorsSetup before using NumeratorsToMasterNumerators.  You must also RunJacobiSetup and set up nBasis either by hand or with FindSomeBasisOfNumerators.";

RunNumeratorsToMasterNumeratorsSetup::usage="RunNumeratorsToMasterNumeratorsSetup sets up NumeratorsToMasterNumerators so that you can use it.  You must RunJacobiSetup and set up nBasis either by hand or with FindSomeBasisOfNumerators before using RunNumeratorsToMasterNumeratorsSetup.";


EvaluateJacobi::usage="EvaluateJacobi[jacobi, RepNumFunc].  Here jacobi is a list of numerators appearing in a Jacobi relation {num[1][...],-num[2][...],-num[7][...]}.  For example jacobi could come from JacobiRelationsOnNumerators.  RepNumFunc[expression] (replace numerator function) is a function that turns every element of nBasis into an actual expression in terms of Mandelstams.  RepNumFunc is probably the function that inserts your actual ansatz.  EvaluateJacobi returns the CoefficientArrays of the Jacobi relation over the MandelstamBasis of the pinched graph associated with the Jacobi relation.  Before calling EvaluateJacobi you need to RunNumeratorsToMasterNumeratorsSetup first.";


RunGraphSymSetup::usage="RunGraphSymSetup sets up the graph symmetries by populating GraphSymRelationsOnNumerators and GraphSymRelationsOnNumeratorsSorted.  NumeratorsToMasterNumerators must have been set up first.";

GraphSymRelationsOnNumerators::usage="GraphSymRelationsOnNumerators is a list of graph symmetry relations on the numerators num[...][...].  Each entry in the list represents an equation of the form num[][]+num[][]==0.  GraphSymRelationsOnNumerators is populated by RunGraphSymSetup.";

GraphSymRelationsOnNumeratorsSorted::usage="GraphSymRelationsOnNumeratorsSorted takes GraphSymRelationsOnNumerators and looks at what they would look like after using NumeratorsToMasterNumerators to go to the basis of numerators in nBasis.  The graph symmetry equations are then sorted by length.  Graph symmetry equations that are trivially satisfied are discarded.  GraphSymRelationsOnNumeratorsSorted returns the the graph symmetries in terms of the original numerators (before NumeratorsToMasterNumerators was applied) so that you can identify what graph kinematics to use when solving the graph symmetry.  Sorting the equations means that you will keep the sparsest ones when finding the linearly independent equations using my solver.  NumeratorsToMasterNumerators must be set up before using GraphSymRelationsOnNumeratorsSorted.  GraphSymRelationsOnNumeratorsSorted is populated by RunGraphSymSetup.";

EvaluateGraphSymRelation::usage="EvaluateGraphSymRelation[GraphSymRelation,RepNumFunc].  Here GraphSymRelation (graph symmetry relation) is one of the entries in GraphSymRelationsOnNumerators or GraphSymRelationsOnNumeratorsSorted.  RepNumFunc[expression] (replace numerator function) is a function that turns every element of nBasis into an actual expression in terms of Mandelstams.  RepNumFunc is probably the function that inserts your actual ansatz.  EvaluateGraphSymRelation returns the CoefficientArrays of the GraphSymRelation over the MandelstamBasis of the graph associated with GraphSymRelation.  Before calling EvaluateGraphSymRelation you need to RunGraphSymSetup first.";


RepListIncludingSigns::usage="RepListIncludingSigns[ListOld,ListNew] generates a set of replacement rules to replace the elements of ListOld with those of ListNew.  You have to be very careful when using RepListIncludingSigns because ListOld can only contain things like p[3] or -p[7] not p[1]-p[2].";


in::usage="in is the Head used by UnorderedCubicGraphs and OrderedCubicGraphs for labeling internal momenta like in[1], in[7]...";


UnorderedCubicGraphs::usage="UnorderedCubicGraphs[n,MomentumHead] returns the n-pt unordered cubic tree graphs.  MomentumHead is the Head (like p or k) used for external momenta.  'in' is the Head used for internal momenta.  UnorderedCubicGraphs is memoized.";

OrderedCubicGraphs::usage="OrderedCubicGraphs[n,MomentumHead] returns the n-pt unordered cubic tree graphs.  MomentumHead is the Head (like p or k) used for external momenta.  'in' is the Head used for internal momenta.  OrderedCubicGraphs is memoized.";


BlowUpOrderedCuts::usage="BlowUpOrderedCuts[VerGraph] returns the set of graphs obtained by blowing up all of the blobs in VerGraph into (color) ordered cubic trees.  Here VerGraph represents what I call a 'graph of cuts' meaning that every drawn edge is cut/on-shell.  So the blobs just correspond to vertices.";

BlowUpUnorderedCuts::usage="BlowUpUnorderedCuts[VerGraph] returns the set of graphs obtained by blowing up all of the blobs in VerGraph into unordered cubic trees.  Here VerGraph represents what I call a 'graph of cuts' meaning that every drawn edge is cut/on-shell.  So the blobs just correspond to vertices.";


GetBlownUpGraphMomConsAndUncutProps::usage="GetBlownUpGraphMomConsAndUncutProps[GraphOfCuts,BlownUpGraph] is used in performing cuts where GraphOfCuts and BlownUpGraph are VerGraphs.  In GraphOfCuts every drawn propagator is cut so blobs are just the vertices.  BlownUpGraph is one of the graphs obtained by blowing up the blobs/vertices in GraphOfCuts using BlowUpOrderedCuts or BlowUpUnorderedCuts.  GetBlownUpGraphMomConsAndUncutProps returns {MomConsRules,UncutMomenta} where UncutMomenta are the momenta inserted in GraphOfCuts to blow it up into BlownUpGraph.  MomConsRules is the set of momentum conservation rules needed to eliminate UncutMomenta and express them in terms of the other momenta in the graph.  GetBlownUpGraphMomConsAndUncutProps will need to be followed up with GetGraphOfCutsKinematics so that CutLHS and CutRHS will be phrased in the same variables.";

GetGraphOfCutsKinematics::usage="GetGraphOfCutsKinematics[GraphOfCuts].  GraphOfCuts is a VerGraph where every drawn propagator is cut so blobs are just the vertices.  GetGraphOfCutsKinematics returns {MomConsRules,dRules} where dRules enfoce that all drawn propagators are on-shell.";


MandelstamBasisOfGraphOfCuts::usage="MandelstamBasisOfGraphOfCuts[GraphOfCuts] returns the basis of momentum dot products associated with the graph GraphOfCuts where all of the drawn propagators (internal and external) have been put on-shell.  This is the basis used by GetGraphOfCutsKinematics.";


GraphOfCutsHasBlownUpZeroPropQ::usage="GraphOfCutsHasBlownUpZeroPropQ[GraphOfCuts] returns True if there is any way of blowing up GraphOfCuts into *unordered* tree graphs that produces a propagator that is identically zero.  GraphOfCutsHasBlownUpZeroPropQ is useful for filtering out unphysical or subtle cuts.  This filters out graphs that blow up to include tadpoles or bubbles on external legs as well as other subtle cuts.";


GenInternalLabelingsOfGraphOfCuts::usage="GenInternalLabelingsOfGraphOfCuts[GraphOfCuts].  GraphOfCuts is a VerGraph representing a cut where every drawn propagator is on-shell and blobs are represented by vertices in GraphOfCuts.  GenInternalLabelingsOfGraphOfCuts produces a list of GraphsOfCuts where every blob goes over (n-3)! color orderings.  GenInternalLabelingsOfGraphOfCuts will thus generate all of the color ordered cuts you need to do for a theory that has tree-level color-kinematics duality like NLSM or YM.";


GenNkMCGraphsOfCuts::usage="GenNkMCGraphsOfCuts[k] generates all of the color ordered, N^k max cuts.  The underlying theory is assumed to have color-kinematics duality at tree level so that you only need to take (n-3)! permutations of the legs on each blob.  The output of GenNkMCGraphsOfCuts is a list of VerGraphs representing GraphsOfCuts where every drawn line is cut (on-shell) and blobs are represented by vertices.  Subtle/non-physical cuts are not produced.  GenNkMCGraphsOfCuts is memoized for speed but the cache is cleared by RunGraphSetup.";


CutLHS::usage="CutLHS[GraphOfCuts, AmpHead] returns the cut constructed from tree graphs.  The cut is represented by GraphOfCuts where every drawn propagator is cut and the vertices represent blobs.  The cut is returned in the basis given by MandelstamBasisOfGraphOfCuts[GraphOfCuts].  AmpHead is the Head used to evaluate the tree amplitudes.  So for exmaple AmpHead could be NLSM so that AmpHead[{p[1],p[2]...}] returns the NLSM tree amplitude.";


CutRHS::usage="CutRHS[GraphOfCuts,NumRepFunc] returns the cut constructed from the numerators of CubicBasisGraphs.  The cut is represented by GraphOfCuts where every drawn propagator is cut and the vertices represent blobs.  The cut is returned in the basis given by MandelstamBasisOfGraphOfCuts[GraphOfCuts].  NumRepFunc is a function that takes any num[...][...] and converts it into the desired representation.  For example NumRepFunc=Identity will return the cut in terms of num[][] which can be useful for debugging.  Alternatively something like NumRepFunc=RepNumAnsatz[NumeratorsToMasterNumerators[#]]& will explicitly evaluate the cut in terms of some ansatz where NumeratorsToMasterNumerators first converts any num[][] into a basis num[][] and RepNumAnsatz replaces the basis num[][] with an actual ansatz.";


Begin["`Private`"]


(* ::Text:: *)
(*TODO:  change all replacement rules (MomConsRepRules, dRepRules, ...) to be functions that act on graphs.  For example change MomConsRepRules[VerGraph,MomToSolveFor:{}] to MomConsRepRules[MomToSolveFor:{}][VerGraph] so that you can call these in a simple way using //.*)


(* ::Subsection:: *)
(*Main code*)


(* ::Subsection:: *)
(*Graph generation and general things*)


GenCubicTreeTopologies[n_,MomHead_]:=Module[{a=MomHead,groupings,f,i,graphs,ExtVerts},
groupings=Groupings[Array[a,n-1],f->{2,Orderless}];
graphs=Table[i=n-1;
group//.f[g_[x_],h_[y_]]:>(i++;Sow[VerF[g[x],h[y],-a[i]]];a[i])//Reap//Last//Last
,{group,groupings}]/.-a[i]:>a[i];

ExtVerts=Join[VerF[-#]&/@Array[a,n-1],{VerF[-a[i]]}];

graphs=VerGraph[Join[#,ExtVerts]]&/@graphs;

DeleteDuplicatesBy[graphs,VerGraphHashCode](*Faster because it only calculates the canonical graph once.  Also it's fine to use the normal CanonicalGraph here.*)
];


MultigraphDisambiguateHelper[MMAGraph_]:=(If[!MultigraphQ[MMAGraph],Return[MMAGraph]];
Module[{edgeList,gatheredEdges,duplicateEdges,inequivalentEdges,i=0,h},
edgeList=EdgeList[MMAGraph];
gatheredEdges=edgeList/.UndirectedEdge[x__]:>UndirectedEdge@@Sort[{x}]//Gather;(*UndirectedEdge isn't Orderless so I Sort it.*)
inequivalentEdges=gatheredEdges//First/@#&//Flatten;
duplicateEdges=gatheredEdges//Rest/@#&//Flatten;
Graph[Flatten[Join[inequivalentEdges,Table[i++;{UndirectedEdge[edge[[1]],h[i]],UndirectedEdge[h[i],edge[[2]]]},{edge,duplicateEdges}]]]]
]);

MultigraphDisambiguate[MMAGraph_]:=FixedPoint[MultigraphDisambiguateHelper,MMAGraph];


(*Clear[VerGraphHashCode];*)

DefineVerGraphHashCode:=(
VerGraphHashCode[VerGraph[VerList_]]:=VerGraphHashCode[VerGraph[VerList]]=VerGraph[VerList]//VerGraphToMMAGraph//MultigraphDisambiguate//CanonicalGraph;
);


GetExtVerts[VerGraph[VerList_]]:=Select[VerList,#/.VerF[x__]:>Length[{x}]==1&];


(*I wouldn't trust this code with 2pt functions.*)
GenCubicTopologies[n_,l_,MomHead_]:=Module[{a=MomHead,CloseLoop,graphs,extLegs,intLegs},
(*Close two external lines to form a loop*)
CloseLoop[graphs_]:=DeleteDuplicatesBy[Flatten[Table[VerGraph[DeleteElements[Identity@@graph,pair]/.-pair[[1,1]]->pair[[2,1]]],{graph,graphs},{pair,Subsets[GetExtVerts[graph],{2}]}]],VerGraphHashCode];

(*Close lines repeatedly*)
graphs=Nest[CloseLoop,GenCubicTreeTopologies[n+2l,a],l];

(*Relabel graphs so external legs are labeled a[1],a[2]...a[n] and internal legs are labeled a[n+1], a[n+2]...*)
Table[
extLegs=gr//GetExtVerts//Identity@@#&/@#&//#/.-x_:>x&;
intLegs=Complement[gr//Identity@@#&//Sequence@@#&/@#&//#/.-x_:>x&//DeleteDuplicates,extLegs];
gr/.MapThread[Rule,{Join[extLegs,intLegs],Array[a,2n-3+3l]}]
,{gr,graphs}]
];


(* ::Input:: *)
(*(*MomConsRepRules[VerGraph[VerList_],MomToSolveFor_:{}]:=Module[{IntMomenta,ExtMomenta},*)
(*ExtMomenta=GetExtVerts[VerGraph[VerList]]/.VerF[x_]:>x/.-x_:>x;*)
(*IntMomenta=Complement[VerList/.-x_:>x/.VerF->List//Flatten//DeleteDuplicates,ExtMomenta];*)
(*Reduce[DeleteElements[VerList,GetExtVerts[VerGraph[VerList]]]/.VerF[x__]:>Total[{x}]==0,Join[IntMomenta,{Last[ExtMomenta]},MomToSolveFor]]//ToRules(*Reduce[eqns,vars] seems to solve eqns in terms of the *reverse* order of vars.  So you want to solve for MomToSolveFor first and then if you can solve for the internal momenta in terms of the external momenta.*)*)
(*];*)*)


(* ::Input:: *)
(*(*MomConsRepRules[VerGraph[VerList_],MomToSolveFor_:{}]:=Module[{IntMomenta,ExtMomenta},*)
(*ExtMomenta=GetExtVerts[VerGraph[VerList]]/.VerF[x_]:>x/.-x_:>x;*)
(*IntMomenta=Complement[VerList/.-x_:>x/.VerF->List//Flatten//DeleteDuplicates,ExtMomenta];*)
(*Reduce[DeleteElements[VerList,GetExtVerts[VerGraph[VerList]]]/.VerF[x__]:>Total[{x}]==0,Join[Complement[IntMomenta,MomToSolveFor],{Last[ExtMomenta]},MomToSolveFor]]//ToRules(*Reduce[eqns,vars] seems to solve eqns in terms of the *reverse* order of vars.  So you want to solve for MomToSolveFor first and then if you can solve for the internal momenta in terms of the external momenta.*)*)
(*];*)*)


MomConsRepRules[VerGraph[VerList_],MomToSolveFor_:{}]:=Module[{IntMomenta,ExtMomenta},
ExtMomenta=GetExtVerts[VerGraph[VerList]]/.VerF[x_]:>x/.-x_:>x;
IntMomenta=Complement[VerList/.-x_:>x/.VerF->List//Flatten//DeleteDuplicates,ExtMomenta];
Reduce[DeleteElements[VerList,GetExtVerts[VerGraph[VerList]]]/.VerF[x__]:>Total[{x}]==0,Join[Complement[Join[IntMomenta,{Last[ExtMomenta]}],MomToSolveFor],MomToSolveFor]]//ToRules(*Reduce[eqns,vars] seems to solve eqns in terms of the *reverse* order of vars.  So you want to solve for MomToSolveFor first and then if you can solve for the internal momenta in terms of the external momenta.*)
];


dRepRules[VerGraph[VerList_],AdditionalConstraints_:{},MomToSolveFor_:{}]:=Module[{ExtMomenta},
ExtMomenta=GetExtVerts[VerGraph[VerList]]/.VerF[x_]:>x/.-x_:>x;
Reduce[Join[(d2[#]==0&/@ExtMomenta),AdditionalConstraints]/.MomConsRepRules[VerGraph[VerList],MomToSolveFor],MomToSolveFor]//ToRules
];


RepGraphKin[graph_][expr_]:=expr/.MomConsRepRules[graph]/.dRepRules[graph];


MandelstamBasisOfGraph[VerGraph[VerList_]]:=Module[{vars},
vars=VerList//#/.VerF->List&//#/.-x_:>x&//Flatten//DeleteDuplicates;
Outer[d,vars,vars]//Flatten//RepGraphKin[VerGraph[VerList]]//#/.d[x__]:>Sow[d[x]]&//Reap//Last//First//DeleteDuplicates
];


TadpoleQ[graph_]:=MomConsRepRules[graph]//Last/@#&//MemberQ[0];


(* ::Text:: *)
(*Tadpole:  momentum conservation sets an internal leg to zero*)
(*Bubble on external leg or tadpole:  on-shell conditions set \ell^2 to zero for some internal leg.*)


TadpoleOrBubbleOnExtLegQ[VerGraph[VerList_]]:=Module[{ExtMomenta,IntMomenta},
ExtMomenta=GetExtVerts[VerGraph[VerList]]/.VerF[x_]:>x/.-x_:>x;
IntMomenta=Complement[VerList/.-x_:>x/.VerF->List//Flatten//DeleteDuplicates,ExtMomenta];
(d2/@IntMomenta)/.MomConsRepRules[VerGraph[VerList]]/.dRepRules[VerGraph[VerList]]/.d[0,_]:>0//Expand//MemberQ[0]
];


(* ::Text:: *)
(*I'd really prefer a function that just replaces all of the d[x,y] for a graph instead of having to first do MomConsRepRules and then dRepRules.  Of course if you have d[e,p] then you'd have to add this in too so maybe this solution will do for now...*)


(* ::Text:: *)
(*Graph generation checks:*)
(*4 legs, 4 loops (about 100 seconds  total)*)
(*5159*)
(*1391 -- no tadpoles*)
(*815 -- no tadpoles or bubbles on external legs*)
(**)
(*3 legs, 3 loops*)
(*147*)
(*45 -- no tadpoles*)
(*17 -- no tadpoles or bubbles on external legs*)


(* ::Input:: *)
(*(*AddLVecHeads[a];*)
(*(grs=GenCubicTopologies[3,3,a];grs//Length)//AbsoluteTiming*)
(*Select[grs,(!TadpoleQ[#])&]//Length//AbsoluteTiming*)
(*Select[grs,(!TadpoleOrBubbleOnExtLegQ[#])&]//Length//AbsoluteTiming*)
(*Clear[grs]*)*)


GraphSignature[VerGraph[VerList_]]:=VerList/.VerF->List//Signature/@#&//Times@@#&;


VerGraphToMMAGraph[VerGraph[VerList_]]:=Graph[UndirectedEdge@@@Table[VerList[[#]]&/@First/@Position[VerList,vert],{vert,VerList/.-x_:>x/.VerF->List//Flatten//DeleteDuplicates}]];(*Needs to capture tadpoles and the 3pt correctly so you have to use this ugly solution*)


VerGraphToEdgeGraph[VerGraph[VerList_]]:=Table[DirectedEdge[VerList[[FirstPosition[VerList,vert,Missing["NotFound"],{2}][[1]]]],VerList[[FirstPosition[VerList,-vert,Missing["NotFound"],{2}][[1]]]]]->vert,{vert,VerList/.-x_:>x/.VerF->List//Flatten//DeleteDuplicates}]//EdgeGraph;(*The funny syntax on FirstPosition is so that you pick up vert and not -vert by mistake.*)


(* ::Input:: *)
(*(*EdgeGraphToUndirectedMMAGraph[EdgeGraph[edgeList_]]:=edgeList//First/@#&//# /.DirectedEdge->UndirectedEdge&//Graph;*)*)


EdgeGraphToVerGraph[EdgeGraph[edgeList_]]:=edgeList//First/@#&//VertexList//VerGraph;


PrintEdgeGraph[EdgeGraph[edgeList_]]:=Module[{extVerts,nn},
extVerts=EdgeGraph[edgeList]//EdgeGraphToVerGraph//GetExtVerts;
nn=Length[extVerts];

GraphPlot[edgeList//Labeled[#[[1]],#[[2]]//ToString]&/@#&,VertexCoordinates->MapThread[Rule,{extVerts,Table[N[{Cos[2Pi i/nn],Sin[2Pi i/nn]}],{i,1,nn}]}]]
];


(* ::Input:: *)
(*(*PrintEdgeGraphOrdered[EdgeGraph[edgeList_],OrderedEdgeLabels_]:=Module[{nn},*)
(*nn=Length[OrderedEdgeLabels];*)
(**)
(*GraphPlot[edgeList//Labeled[#[[1]],#[[2]]//ToString]&/@#&,VertexCoordinates->MapThread[Rule,{OrderedEdgeLabels//VerF[-#]&/@#&,Table[N[{Cos[2Pi i/nn],Sin[2Pi i/nn]}],{i,1,nn}]}]]*)
(*];*)*)


(* ::Subsection:: *)
(*Evil test graphs*)


(* ::Input:: *)
(*(*VertexShorthand[expr_]:=expr/.v[x__]:>VerF@@a/@{x}/.a[x_]:>-a[Abs[x]]/;x<0;*)
(**)
(*tadpole=VerGraph[{v[1,2,-1],v[-2]}]//VertexShorthand;*)
(*tadpole2=VerGraph[{v[1,2,3,-2,-1],v[-3]}]//VertexShorthand;*)
(*singlevertex=VerGraph[{v[1,2,3],v[-1],v[-2],v[-3]}]//VertexShorthand;*)
(*evilbubble=VerGraph[{v[1,2,3],v[-2,4,-5],v[-4,6,5],v[-6,9,10],v[7,-9,8],v[-3,-7,-8],v[-1],v[-10]}]//VertexShorthand;*)
(*evilbubble2=VerGraph[{v[1,3,-4,5],v[-3,4,-5,2],v[-1],v[-2]}]//VertexShorthand;*)
(*evilbubble3=VerGraph[{v[1,3,-4,5],v[-3,4,-5,6],v[-6,7,8,-9],v[9,-8,-7,10],v[-1],v[-10]}]//VertexShorthand;*)
(*testtree=GenCubicTopologies[8,0,a][[2]];*)
(*singlebubble=VerGraph[{v[-1],v[1,2,3],v[-2,-3,4],v[-4]}]//VertexShorthand;*)*)


(* ::Input:: *)
(*(*tadpole//VerGraphToMMAGraph//MultigraphQ ;(*False*)*)
(*tadpole2//VerGraphToMMAGraph//MultigraphQ ;(*True*)*)*)


(* ::Input:: *)
(*(*{tadpole,tadpole2,evilbubble,evilbubble2,evilbubble3}//VerGraphToMMAGraph/@#&;*)*)


(* ::Input:: *)
(*(*BadGraphOfCuts=VerGraph[{VerF[p[1],p[2],-p[7]],VerF[-p[9],-p[5],-p[6],p[7]],VerF[p[9],p[5],-p[10]],VerF[p[10],p[6],-p[11]],VerF[p[11],p[3],p[4]],VerF[-p[1]],VerF[-p[2]],VerF[-p[3]],VerF[-p[4]]}];*)
(**)
(*BadGraphOfCuts//VerGraphHashCode;*)*)


(* ::Subsection:: *)
(*Things for 4pt 3loop N=4*)


(* ::Input:: *)
(*(*BubbleOrTriangleQ[graph_]:=graph//VerGraphToMMAGraph//FindCycle[#,Infinity,All]&//Length/@#&//Min//#<4&;*)
(**)
(*AddLVecHeads[p];*)
(**)
(*Select[GenCubicTopologies[4,3,p],(!TadpoleOrBubbleOnExtLegQ[#])&&(!BubbleOrTriangleQ[#])&]//Length//AbsoluteTiming*)*)


(* ::Subsection:: *)
(*Graph iso code including multigraphs but not tadpoles*)


(* ::Text:: *)
(*Graph iso code won't work well with tadpoles!!!*)


FindVerGraphMorphismsNoMultigraphs[VerGraph[VerListSource_],VerGraph[VerListTarget_],isofunction_]:=Module[{EdGrSource,EdGrTarget,MMAGrSource,MMAGrTarget,GrIsos},

EdGrSource=VerGraphToEdgeGraph[VerGraph[VerListSource]]/.EdgeGraph->Identity;(*(VerS -> VerS) -> eS where S means source*)
EdGrTarget=VerGraphToEdgeGraph[VerGraph[VerListTarget]]/.EdgeGraph->Identity;(*(VerT -> VerT) -> eT where T means target*)

MMAGrSource=EdGrSource//First/@#&//UndirectedGraph;
MMAGrTarget=EdGrTarget//First/@#&//UndirectedGraph;

GrIsos=Sequence[MMAGrSource,MMAGrTarget]//isofunction;(*VerS -> VerT*)

Table[EdGrTarget/.(EdGrSource/.iso)/.(EdGrSource/.iso/.(DirectedEdge[a_,b_]->c_):>(DirectedEdge[b,a]->-c))/.(-x_->y_):>(x->-y),{iso,GrIsos}]
];


DotThisEdge[edgelist_,edge_]:=Module[{vPair,vLeft,vRight,b},
vPair=edgelist//Reverse/@#&//edge/.#&//List@@#&;
vLeft=vPair//First;
vRight=vPair//Last;
edgelist/.(DirectedEdge[vLeft,vRight]->edge)->{DirectedEdge[vLeft,VerF[-edge,b]]->edge,DirectedEdge[VerF[-edge,b],vRight]->b}/.vRight->(vRight/.-edge->-b)//Flatten
];


VerGraphMultigraphDisambiguate[VerGraph[VerList_]]:=Module[{EdList,EdgesToDot},
EdList=VerGraph[VerList]//VerGraphToEdgeGraph//Identity@@#&;
EdgesToDot=EdList//# /.DirectedEdge[x_,y_]:>DirectedEdge@@Sort[{x,y}]&//GatherBy[#,First]&//Select[#,Length[#]>1&]&//Flatten//Last/@#&;
Fold[DotThisEdge,EdList,EdgesToDot]//EdgeGraph//EdgeGraphToVerGraph
];


(* ::Text:: *)
(*TODO:  Just like you had to used FixedPoint in MultigraphDisambiguate, you might need to do it here in VerGraphMultigraphDisambiguate for iterated tadpoles like my tadpole2 example graph.  Also, if you have to do FixedPoint in VerGraphMultigraphDisambiguate then you probably have to be careful to do things right in FindVerGraphIsomorphisms to get rid of all of the extra dots.  You could have two dots on one line.*)


FindVerGraphMorphisms[VerGraph[VerListSource_],VerGraph[VerListTarget_],isofunction_]:=Module[{gr1,gr2,gr1UnDotRules,gr2UnDotRules},
gr1=VerGraph[VerListSource]//VerGraphMultigraphDisambiguate;
gr2=VerGraph[VerListTarget]//VerGraphMultigraphDisambiguate;

If[VerListSource===VerListTarget,gr2=gr1];(*You need this when you're using FindMinimalSpanningAutomorphismsNoMultigraphs to look at multigraphs otherwise the new private labels x$111 from VerGraphMultigraphDisambiguate won't align properly.*)

gr1UnDotRules=gr1//Identity@@#&//#/.VerF[-x_,y_]:>Sow[y->x]&//Reap//Last//Flatten;
gr2UnDotRules=gr2//Identity@@#&//#/.VerF[-x_,y_]:>Sow[y->x]&//Reap//Last//Flatten;

Map[DeleteDuplicates[#/.gr1UnDotRules/.gr2UnDotRules]&,FindVerGraphMorphismsNoMultigraphs[gr1,gr2,isofunction]]
];


(*This might not find the minimal set of automorphism group generators but it should hopefully be complete and I put some effort into cutting down the number of generators.*)
FindMinimalSpanningAutomorphismsNoMultigraphs[MMAGraph_]:=Module[{MMACanonicalGraph,MMACanGrToMMAGrIso,MMAGrToMMACanGrIso,AutoGroup,CyclesToKeep},

MMACanonicalGraph=MMAGraph//CanonicalGraph;
MMACanGrToMMAGrIso=FindGraphIsomorphism[MMACanonicalGraph,MMAGraph]//First;
MMAGrToMMACanGrIso=FindGraphIsomorphism[MMAGraph,MMACanonicalGraph]//First;
AutoGroup=MMACanonicalGraph//GraphAutomorphismGroup;
(*Normally the rank of GraphAutomorphismGroup gives the number of symmetry factors but here I will delete redundant isomorphisms so this won't be true any more.*)
CyclesToKeep=AutoGroup//Identity@@#&;

(*Delete a Cycle (group generator) if it doesn't produce unique isomorphisms (that aren't produced by other Cycles).  I think you can do this by looking at the GroupOrbits.*)
Table[
If[GroupOrbits[CyclesToKeep//DeleteCases[cycle]//PermutationGroup,MMACanonicalGraph//VertexList]===GroupOrbits[AutoGroup,MMACanonicalGraph//VertexList],
CyclesToKeep=CyclesToKeep//DeleteCases[cycle]]
,{cycle,AutoGroup//Identity@@#&}];

Table[MMAGrToMMACanGrIso//PermutationReplace[#,cycle]&/@#&//MMACanGrToMMAGrIso/@#&,{cycle,CyclesToKeep}];

Table[MMAGrToMMACanGrIso//PermutationReplace[#,cycle]&/@#&//MMACanGrToMMAGrIso/@#&,{cycle,CyclesToKeep}]
];


FindVerGraphIsomorphisms[g1_,g2_] := FindVerGraphMorphisms[g1,g2,FindGraphIsomorphism[#1,#2,All]&];

FindVerGraphSpanningAutomorphisms[graph_]:=FindVerGraphMorphisms[graph,graph,FindMinimalSpanningAutomorphismsNoMultigraphs[#1]&];


(* ::Subsection:: *)
(*Tree-level symmetry factor checks*)


(* ::Input:: *)
(*(*AddLVecHeads[a]*)*)


(* ::Input:: *)
(*(*Table[{GenCubicTreeTopologies[n,a]//FindVerGraphIsomorphisms[#,#]&/@#&//Length/@#&//n!/#&/@#&//Total,(2n-5)!!},{n,4,10}]*)*)


(* ::Subsection:: *)
(*Setup basis graphs*)


RunGraphSetup[n_,l_,MomentumHead_]:=Module[{a=MomentumHead},

ClearAll[GenNkMCGraphsOfCuts,VerGraphHashCode];(*Clear the caches of these functions because if you change CubicBasisGraphs you probably don't need to store the VerGraphHaseCodes any more and you definitely don't want to store GenNkMCGraphsOfCuts.*)

DefineGenNkMCGraphsOfCuts;
DefineVerGraphHashCode;

AddLVecHeads[a];
Print["Added ", a," to LVecHeads"];
CubicBasisGraphs=Select[GenCubicTopologies[n,l,a],(!TadpoleOrBubbleOnExtLegQ[#])&];

CanonicalNumeratorVariables=CubicBasisGraphs//First//Identity@@#&//List@@#&/@#&//Flatten//#/.-x_:>x&//DeleteDuplicates//Sort;

BasisGraphNumeratorPairs=Table[{CubicBasisGraphs[[i]],num[i][CanonicalNumeratorVariables]},{i,1,Length[CubicBasisGraphs]}];

BasisMMAGraphNumeratorPairs=Table[{elm[[1]]//VerGraphHashCode,elm[[2]]/.num[x_][y_]:>num[x]},{elm,BasisGraphNumeratorPairs}];
];


(* ::Subsection:: *)
(*Generate functional Jacobis*)


PinchEdgeHelper[VerGraph[VerList_],EdGrEdgeToPinch_]:=Module[{VerL,VerR,NewVer,SharedEdge},
SharedEdge=EdGrEdgeToPinch[[2]];
VerL=EdGrEdgeToPinch//List@@#[[1,1]]&;
VerR=EdGrEdgeToPinch//List@@#[[1,2]]&;
NewVer=Join[VerL/.{a___,SharedEdge,b___}:>{Reverse[{b}],a},VerR/.{a___,-SharedEdge,b___}:>{b,Reverse[{a}]}]//Flatten//VerF@@#&;
DeleteCases[VerGraph[VerList],VerF@@VerL,Infinity]/.VerF@@VerR->NewVer
];


(*Here the edge must be a *positive* momentum label like p[4] not -p[3]. <-- This shouldn't be true anymore.*)
PinchSpecifiedEdge[VerGraph[VerList_],edge_]:=Module[{EdGraph,EdGrEdgeToPinch,ed},
ed=edge/.-x_:>x;
EdGraph=VerGraphToEdgeGraph[VerGraph[VerList]];EdGrEdgeToPinch=EdGraph//Identity@@#&//Select[#,Last[#]===edge&]&//First;PinchEdgeHelper[VerGraph[VerList],EdGrEdgeToPinch]
];


PinchOneEdge[VerGraph[VerList_]]:=Module[{EdGraph,EdList,SharedEdge,VerL,VerR,NewVer,ExtVertList,EdListNoExtVerts,EdListNoExtVertsOrSelfLoops},
EdGraph=VerGraph[VerList]//VerGraphToEdgeGraph;
EdList=EdGraph//Identity@@#&;
ExtVertList=GetExtVerts[VerGraph[VerList]]//Identity@@#&/@#&//#/.-x_:>x&;
EdListNoExtVerts=Select[EdList,!MemberQ[ExtVertList,Last[#]]&];
EdListNoExtVertsOrSelfLoops=Select[EdListNoExtVerts,(!MatchQ[#,DirectedEdge[VerF[x__],VerF[x__]]->y_])&];

EdListNoExtVertsOrSelfLoops//PinchEdgeHelper[VerGraph[VerList],#]&/@#&//DeleteDuplicatesBy[#,VerGraphHashCode]&
];


PinchOneEdgeOnSetOfGraphs[VerGraphList_]:=VerGraphList//PinchOneEdge/@#&//Flatten//DeleteDuplicatesBy[#,VerGraphHashCode[#]&]&;


BlowUpVertexIntoJacobiTriplet[VerGraph[VerList_]]:=Module[{ax},
ax=EdgeInsertedInJacobi;
{VerGraph[VerList]/.VerF[a1_,a2_,a3_,a4_]:>Sequence[VerF[a1,a2,ax],VerF[-ax,a3,a4]],VerGraph[VerList]/.VerF[a1_,a2_,a3_,a4_]:>Sequence[VerF[a1,a3,ax],VerF[-ax,a4,a2]],VerGraph[VerList]/.VerF[a1_,a2_,a3_,a4_]:>Sequence[VerF[a1,a4,ax],VerF[-ax,a2,a3]]}
];(*Be careful about the labels on how you blow up the vertex.  Later on we will use EdgeInsertedInJacobi to easily identify the inserted leg to eliminate it with momentum conservation.*)


FindVerGraphNumerator[VerGraph[InputVerList_]]:=Module[{CanMMAGraph,ni,basisVerGraph,iso},
CanMMAGraph=VerGraph[InputVerList]//VerGraphHashCode;
ni=CanMMAGraph/.(BasisMMAGraphNumeratorPairs//Rule@@#&/@#&);
basisVerGraph=ni/.(BasisGraphNumeratorPairs//Reverse/@#&//#/.num[x_][y_]:>num[x]&//Rule@@#&/@#&);
(*Find isomorphism from basis graph to the input graph*)
iso=FindVerGraphIsomorphisms[basisVerGraph,VerGraph[InputVerList]]//First;
GraphSignature[basisVerGraph/.iso]GraphSignature[VerGraph[InputVerList]] (ni[CanonicalNumeratorVariables]/.iso)(*Make sure the signs make sense.  Isomorphism phi maps from source graph (g_sr) to target graph (g_ta) and sgn(g) means the signature of graph g.  g_ta = sgn(g_ta) sgn(phi(g_sr)) phi(g_sr)*)
];


RunJacobiSetup[IncludeTadpoleBELJacobis_]:=(
EdgeInsertedInJacobi=(CanonicalNumeratorVariables//First//Head)[CanonicalNumeratorVariables//Length//#+1&];

JacobiRelationsOnGraphs=CubicBasisGraphs//PinchOneEdgeOnSetOfGraphs//BlowUpVertexIntoJacobiTriplet/@#&//Select[#,(!TadpoleOrBubbleOnExtLegQ[#])&]&/@#&;

If[IncludeTadpoleBELJacobis==False,JacobiRelationsOnGraphs=JacobiRelationsOnGraphs//Select[#,Length[#]==3&]&;];

JacobiRelationsOnNumerators=Module[{EliminateInsertedLeg},
Table[
EliminateInsertedLeg=EdgeInsertedInJacobi->(graph//Identity@@#&//SelectFirst[#,(!FreeQ[#,-EdgeInsertedInJacobi]&)]&//List@@#&//Total//#+EdgeInsertedInJacobi&);graph//FindVerGraphNumerator//#/.EliminateInsertedLeg&
,{jacobi,JacobiRelationsOnGraphs},{graph,jacobi}]];(*Find the Jacobi relations between the numerators and eliminate the inserted leg with momentum conservation.*)

SimplifiedJacobis=JacobiRelationsOnNumerators/.num[x_][y_]:>num[x]/.a_*num[x_]:>num[x]//Sort/@#&//DeleteCases[#,{x_,x_,x_}]&//DeleteCases[#,{x_,x_}]&;(*Want to ignore self Jacobis that don't help determine a basis*)
);


(* ::Subsection:: *)
(*Find a (possibly overcomplete) basis of numerators and solve for all numerators in terms of basis*)


FindConstructibleNumerators[nums_]:=Module[{CurrentNums},
CurrentNums=nums;
Table[If[Length[#]==1,AppendTo[CurrentNums,#[[1]]];]&@Fold[DeleteCases,jac,CurrentNums],{jac,SimplifiedJacobis}];(*If you know all but one entry in a Jacobi then you can use the Jacobi to determine the last entry.*)
CurrentNums//Sort
];


FindSomeBasisOfNumerators:=(
nBasis=Module[{AllNumerators,NumBasis},
AllNumerators=Array[num,BasisGraphNumeratorPairs//Length];NumBasis=AllNumerators;

Table[If[FixedPoint[FindConstructibleNumerators,DeleteCases[NumBasis,NumToDelete]]===AllNumerators,NumBasis=DeleteCases[NumBasis,NumToDelete]];
,{NumToDelete,(*AllNumerators*)(*Reverse[AllNumerators]*)(*RandomSample[AllNumerators]*)Reverse[AllNumerators]}];
NumBasis
];
);
(*At 4pt 3loop something interesting happens.  If you choose RandomSample of AllNumerators you will sometimes get a basis of 2 graphs and sometimes a basis of 3 graphs.  In particular, you can generate all of the numerators from the tennis court and another specific graph that has a triangle.  The order in which you try to elliminate the numerators matters!  This algorithm may generate an overcomplete basis.  Tests:  at tree level the basis should just be the half ladder and at 1 loop it should just be the n-gon.*)


SimplifiedNumsToVerGraphs[NumList_]:=NumList/.(BasisGraphNumeratorPairs//#/.num[x_][_]:>num[x]&//Reverse/@#&//Rule@@#&/@#&);


RepListIncludingSigns[ListOld_,ListNew_]:=Module[{signs},
signs=ListOld//Sign/@#&//#/.Sign[_]:>1&;
RepList[signs*ListOld,signs*ListNew]
];(*You have to be very careful in using this because you can't have something like p[1]-p[2] in ListOld only stuff like p[1] or -p[2].*)


NumeratorsToMasterNumerators[expr_]:=Fold[#1/.#2[[1]][list_]:>(#2[[2]]/.RepList[CanonicalNumeratorVariables,list])&,expr,NumJacobiRepTuples];


RunNumeratorsToMasterNumeratorsSetup:=(
NumJacobiRepTuples={};(*{{num[2],-num[4][{...}]+num[4][{...}]},...}*)

Module[{EnumeratableNumerators,JacobisWithUnknownNums,SolvableJacobiIndex,SolvableJacobiIndices,SolvableJacobiLengths,NewNum,nSolved,nRest,CurrentLabeling,nRestRelabeled,nSolvedSign,IndexToDrop},
EnumeratableNumerators=nBasis;
Do[
JacobisWithUnknownNums=JacobiRelationsOnNumerators/.num[x_][y_]:>num[x]/.a_*num[x_]:>num[x]// Fold[DeleteCases,#,EnumeratableNumerators]&/@#&;

(*Old naive method*)
(*SolvableJacobiIndex=JacobisWithUnknownNums//Length/@#&//FirstPosition[1]//Identity@@#&;*)
(*Find the first Jacobi where you know all but one entry*)

(*Find the Jacobis where you know all but one entry.  These Jacobis are solvable.  Next sort them by whichever will produce the shortest numerator and solve that one.*)
SolvableJacobiIndices=JacobisWithUnknownNums//Length/@#&//Position[1]//Flatten;
SolvableJacobiLengths=JacobiRelationsOnNumerators[[SolvableJacobiIndices]]//NumeratorsToMasterNumerators//#/.Plus->List&//Flatten/@#&//#/.num[x_][y_]:>num[x]&//#/.a_*num[x_]:>num[x]&//Length/@#&;
SolvableJacobiIndex=MapThread[List,{SolvableJacobiIndices,SolvableJacobiLengths}]//SortBy[Last]//First//First;

NewNum=JacobisWithUnknownNums[[SolvableJacobiIndex]]//First;

nSolved=Select[JacobiRelationsOnNumerators[[SolvableJacobiIndex]],!FreeQ[#,NewNum]&]//Identity@@#&;
nSolvedSign=Sign[nSolved/.num[_][_]:>1];
nSolved=nSolvedSign*nSolved;(*Flip the sign of nSolved if needed*)

nRest=-nSolvedSign*Total[Select[JacobiRelationsOnNumerators[[SolvableJacobiIndex]],FreeQ[#,NewNum]&]];(*Be carefule about the sign*)

CurrentLabeling=nSolved/.num[x_][y_]:>y;
(*CurrentLabeling will involve one entry that looks schematically like -p1+p2 because you used momentum conservation to eliminate EdgeInsertedInJacobi.  You need to map CurrentLabeling to CanonicalNumeratorVariables but omit this funny -p1+p2 entry.*)
IndexToDrop=CurrentLabeling//FirstPosition[#,_?(Head[#]===Plus&)]&;

nRestRelabeled=nRest/.RepListIncludingSigns[Drop[CurrentLabeling,IndexToDrop],Drop[CanonicalNumeratorVariables,IndexToDrop]]//NumeratorsToMasterNumerators;

AppendTo[NumJacobiRepTuples,{NewNum,nRestRelabeled}];

AppendTo[EnumeratableNumerators,NewNum];
,Length[CubicBasisGraphs]-Length[nBasis]](*Do this enough times until you can rewrite every numerator in terms of the basis*)
];
);


(* ::Text:: *)
(*For solving graph sym relations...*)
(*You're trying to enforce graph sym relations on graph G.  Take the kinematics from graph G and enforce them on the numerators that come from your basis (G1, G2...).*)
(**)
(*For solving Jacobi relations...*)
(*Take a single pinch to produce graph G.  Take the kinematics from graph G and apply it to all of the blown up graphs that constitute the Jacobi triplet that you're trying to solve.  I think you'll still have some kinematics left over to specify like you'll still have some momentum conservation to do.*)


(* ::Subsection:: *)
(*Create a way to actually evaluate a Jacobi relation on numerators*)


EvaluateJacobi[jac_,RepNumFunc_]:=Module[{PinchedGraphLabeling,index,InitLabeling,PinchedGraph},
InitLabeling=jac//First//#/.num[_][y_]:>Sow[y]&//Reap//Last//First//First;index=InitLabeling//FirstPosition[#,_?(Head[#]===Plus&)]&//First;
PinchedGraphLabeling=InitLabeling//#/.InitLabeling[[index]]->EdgeInsertedInJacobi&;
PinchedGraph=CubicBasisGraphs[[jac//First//#/.num[x_][_]:>Sow[x]&//Reap//Last//First//First]]/.RepList[CanonicalNumeratorVariables,PinchedGraphLabeling]//PinchSpecifiedEdge[#,EdgeInsertedInJacobi]&;
jac//Total//NumeratorsToMasterNumerators//RepNumFunc//RepGraphKin[PinchedGraph]//CoefficientArrays[#,MandelstamBasisOfGraph[PinchedGraph]]&
];


(* ::Subsection:: *)
(*Functional relations for numerators respecting graph syms*)


RunGraphSymSetup:=(
GraphSymRelationsOnNumerators=Module[{gr,numer},
Table[
gr=elm[[1]];
numer=elm[[2]];
Table[numer-GraphSignature[gr]GraphSignature[gr/.iso](numer/.iso),{iso,FindVerGraphSpanningAutomorphisms[gr]}]
,{elm,BasisGraphNumeratorPairs}]//Flatten(*Be careful about the signs on GraphSignature*)
];

GraphSymRelationsOnNumeratorsSorted=GraphSymRelationsOnNumerators//{#,#//NumeratorsToMasterNumerators//Length}&/@#&//Select[#,Last[#]>0&]&//SortBy[Last]//First/@#&;

);


EvaluateGraphSymRelation[GraphSymRelation_,RepNumFunc_]:=Module[{graph},
graph=CubicBasisGraphs[[GraphSymRelation//First//#/.num[x_][_]:>Sow[x]&//Reap//Last//First//First]];
GraphSymRelation//NumeratorsToMasterNumerators//RepNumFunc//RepGraphKin[graph]//CoefficientArrays[#,MandelstamBasisOfGraph[graph]]&
];


(* ::Subsection:: *)
(*Unordered and color ordered tree graph generation*)


(*Clear[UnorderedCubicGraphs];*)

UnorderedCubicGraphs[3,MomentumHead_]:=UnorderedCubicGraphs[3,MomentumHead]=GenCubicTreeTopologies[3,MomentumHead];

UnorderedCubicGraphs[n_,MomentumHead_]:=UnorderedCubicGraphs[n,MomentumHead]=Module[{a=MomentumHead},Table[graph/.VerF[x___,ed,y___]:>Sequence[VerF[-a[n]],VerF[ed,in[n-3],a[n]],VerF[x,-in[n-3],y]],{graph,UnorderedCubicGraphs[n-1,MomentumHead]},{ed,Join[Array[a,n-1],Array[in,n-4]//Flatten]}]//Flatten
];


(* ::Input:: *)
(*(*Table[(2n-5)!!,{n,3,8}]*)
(*Table[UnorderedCubicGraphs[n]//Length,{n,3,8}]//AbsoluteTiming*)*)


(* ::Text:: *)
(*So now you need to relabel the graphs to get whatever you want.*)
(*You will have to keep track of internal labels.*)


(*Clear[OrderedCubicGraphs];*)

OrderedCubicGraphs[3,MomentumHead_]:=OrderedCubicGraphs[3,MomentumHead]=GenCubicTreeTopologies[3,MomentumHead];

OrderedCubicGraphs[n_,MomentumHead_]:=OrderedCubicGraphs[n,MomentumHead]=Module[{GoodEdges,edge,vertex,a=MomentumHead},
Table[
edge=a[1];
vertex=graph//Identity@@#&//SelectFirst[#,MatchQ[VerF[___,edge,___]]]&;

GoodEdges=While[UnsameQ[edge , -a[n-1]],
Sow[edge];
vertex=graph//Identity@@#&//SelectFirst[#,MatchQ[VerF[___,edge,___]]]&;
edge=-vertex[[vertex//List@@#&//FirstPosition[edge]//First//Mod[#-1,3(*Cubic vertex*),1]&]];
]//Reap//Last//Flatten;

Join[Table[graph/.VerF[x___,ed,y___]:>Sequence[VerF[-a[n]],VerF[ed,in[n-3],a[n]],VerF[x,-in[n-3],y]],{ed,GoodEdges}],
{graph/.VerF[x___,a[n-1],y___]:>Sequence[VerF[-a[n]],VerF[a[n-1],a[n],in[n-3]],VerF[x,-in[n-3],y]]}
](*The last leg is  weird*)
,{graph,OrderedCubicGraphs[n-1,MomentumHead]}]//Flatten
];


(* ::Input:: *)
(*(*Table[2^(n-2)(2n-5)!!/(n-1)!,{n,3,10}]*)
(*Table[OrderedCubicGraphs[n]//Length,{n,3,10}]//AbsoluteTiming*)*)


(* ::Subsection:: *)
(*MMC code*)


(* ::Input:: *)
(*(*VerGraphSimplifyLabels[VerGraph[VerList_]]:=Module[{SortedEdgeList},*)
(*SortedEdgeList=VerList/.VerF[x__]:>{x}//#/.-x_:>x&//Flatten//DeleteDuplicates//Sort;*)
(*VerGraph[VerList]/.RepList[SortedEdgeList,Array[SortedEdgeList//First//Head,SortedEdgeList//Length]]*)
(*];*)*)


(* ::Text:: *)
(*For a graph of cuts, sort the vertices into two sets*)
(*The first set (the unchanged/prefactor vertices) have length 3 or less.*)
(*The blow up set has length more than 3.*)
(*Loop/Table over the blow up set.  Turn each one into a set of ordered trees and increase the internal leg counter accordingly (i -> i + vertex length -3 where i starts at length of legs in graph).*)
(*Now take outer product by appending the prefactor and then outer producting over every blow up of every graph.*)
(**)
(*So you need to keep track of the uncut propagators (in[1]...) and use momentum conservation to solve for them in terms of the momenta of the GraphOfCuts.  So basically BlowUpOrderedCuts has to keep track of {graph, props, mom cons} for each graph in OrderedCubicGraphs.*)


BlowUpCuts[VerGraph[GraphOfCutsVerList_],TreeGraphGeneratingFunction_]:=Module[{joinGraphs,PrefactorVertices,VerticesToBlowUp,AlmostEdgeList,NodeHead,index,BlownUpGraphs,vList,n,ret,a},
a=GraphOfCutsVerList//First//First//#/.-x_:>x&//Head;(*So 'a' will be the Head of the momentum vectors like 'k' or 'p'.*)
PrefactorVertices=Select[GraphOfCutsVerList,Length[List@@#]<=3&];
VerticesToBlowUp=Complement[GraphOfCutsVerList,PrefactorVertices];
If[VerticesToBlowUp==={},Return[{VerGraph[GraphOfCutsVerList]}]];
AlmostEdgeList=GraphOfCutsVerList//#/.VerF[x__]:>{x}&//Flatten//#/.-x_:>x&;
NodeHead=AlmostEdgeList//First//Head;
index=AlmostEdgeList//#/.NodeHead[x_]:>x&//Max//#+1&;
BlownUpGraphs=Table[
vList=List@@ver;
n=Length[vList];
ret=TreeGraphGeneratingFunction[n,a]/.VerF[-a[x_]]:>(##&[])/.RepList[Join[Array[a,n],Array[in,n-3]],Join[vList,Array[a,n-3,index]]](*/.VerGraph[x_]:>x*);
index+=n-3;
ret
,{ver,VerticesToBlowUp}];
Outer[joinGraphs,##]&@@BlownUpGraphs//Flatten//#/.VerGraph->Identity&//#/.joinGraphs[x__]:>VerGraph[Join[PrefactorVertices,x]]&
];


BlowUpOrderedCuts[VerGraph[GraphOfCutsVerList_]]:=BlowUpCuts[VerGraph[GraphOfCutsVerList],OrderedCubicGraphs];

BlowUpUnorderedCuts[VerGraph[GraphOfCutsVerList_]]:=BlowUpCuts[VerGraph[GraphOfCutsVerList],UnorderedCubicGraphs];


(* ::Text:: *)
(*graphOfCuts, blownUpGraph ---- output the momentum conservation rules to phrase everything in terms of the momenta of the edges of the graphOfCuts and return the uncut propagators*)
(**)
(*All propagators are either internal or external.  All external props need to square to zero.  All internal props are either cut or uncut.  The cut props square to zero.  Need to solve for the uncut props in terms of the cut props.*)


GetBlownUpGraphMomConsAndUncutProps[VerGraph[GraphOfCutsVerList_],VerGraph[BlownUpGraphVerList_]]:=Module[{ExtMomenta,CutMomenta,UncutMomenta,MomConsRules,dRules},
ExtMomenta=GetExtVerts[VerGraph[GraphOfCutsVerList]]/.VerF[x_]:>x/.-x_:>x;
CutMomenta=Complement[GraphOfCutsVerList/.-x_:>x/.VerF->List//Flatten//DeleteDuplicates,ExtMomenta];
UncutMomenta=Complement[BlownUpGraphVerList/.-x_:>x/.VerF->List//Flatten//DeleteDuplicates,CutMomenta,ExtMomenta];

MomConsRules=MomConsRepRules[VerGraph[BlownUpGraphVerList],UncutMomenta];

{MomConsRules,UncutMomenta}
];


GetGraphOfCutsKinematics[VerGraph[GraphOfCutsVerList_]]:=Module[{ExtMomenta,CutMomenta,MomConsRules,dRules},
ExtMomenta=GetExtVerts[VerGraph[GraphOfCutsVerList]]/.VerF[x_]:>x/.-x_:>x;
CutMomenta=Complement[GraphOfCutsVerList/.-x_:>x/.VerF->List//Flatten//DeleteDuplicates,ExtMomenta];
MomConsRules=MomConsRepRules[VerGraph[GraphOfCutsVerList]];
dRules=dRepRules[VerGraph[GraphOfCutsVerList],d2[#]==0&/@CutMomenta];

{MomConsRules,dRules}
];


MandelstamBasisOfGraphOfCuts[VerGraph[GraphOfCutsVerList_]]:=Module[{vars,MomCons,dRules},
{MomCons,dRules}=GetGraphOfCutsKinematics[VerGraph[GraphOfCutsVerList]];vars=GraphOfCutsVerList//#/.VerF->List&//#/.-x_:>x&//Flatten//DeleteDuplicates;Outer[d,vars,vars]//Flatten//#/.MomCons&//#/.dRules&//#/.d[x__]:>Sow[d[x]]&//Reap//Last//First//DeleteDuplicates
];


GraphOfCutsHasBlownUpZeroPropQ[VerGraph[GraphOfCutsVerList_]]:=Module[{BlownUpGraphs,MomConsUncutMom,UncutMomenta,denom,MomCons,dRules},
BlownUpGraphs=BlowUpUnorderedCuts[VerGraph[GraphOfCutsVerList]];(*Blow up unordered cuts to try to see if there are any silly graphs with zero propagators.  Blowing up to ordered cuts won't capture all of the silly graphs.*)
Table[
{MomConsUncutMom,UncutMomenta}=GetBlownUpGraphMomConsAndUncutProps[VerGraph[GraphOfCutsVerList],graph];
{MomCons,dRules}=GetGraphOfCutsKinematics[VerGraph[GraphOfCutsVerList]];
denom=d2/@UncutMomenta//Times@@#&;
denom/.MomConsUncutMom/.MomCons/.dRules//Expand//#/.d[0,_]:>0&
,{graph,BlownUpGraphs}]//MemberQ[0]
];


(* ::Text:: *)
(*For the product of tree amps:*)
(*Turn every graph into a tree amp*)
(*You'll have to normalize the tree amps so that they have the right factors when they factorize on their poles*)
(**)
(*OrderedCubicGraphs[{external leg list}, {internal leg list}]*)


GenInternalLabelingsOfGraphOfCuts[VerGraph[GraphOfCutsVerList_]]:=Module[{PrefactorVertices,VerticesToBlowUp,VerticesWithPermutedLabels,edList,joinGraphs},
PrefactorVertices=Select[GraphOfCutsVerList,Length[List@@#]<=4&];
VerticesToBlowUp=Complement[GraphOfCutsVerList,PrefactorVertices];
If[VerticesToBlowUp==={},Return[{VerGraph[GraphOfCutsVerList]}]];
VerticesWithPermutedLabels=Table[
edList=ver//List@@#&;
edList[[4;;]]//Permutations//Join[edList[[;;3]],#]&/@#&//VerF@@#&/@#&
,{ver,VerticesToBlowUp}];
Outer[joinGraphs,##]&@@VerticesWithPermutedLabels//Flatten//#/.joinGraphs[x__]:>VerGraph[Join[PrefactorVertices,{x}]]&
];(*For each blob in the GraphOfCuts this will do (n-3)! relabelings of the legs.  This should be enough cuts for any theory that has tree-level color-kinematics duality.

This is slight overkill because you only need to permute blobs with more than 4 *internal* legs.  Permuting external legs on a blob won't do anything because you already enforced graph symmetries.*)


(*Clear[GenNkMCGraphsOfCuts];*)

DefineGenNkMCGraphsOfCuts:=(
GenNkMCGraphsOfCuts[k_]:=
GenNkMCGraphsOfCuts[k]=CubicBasisGraphs//Nest[PinchOneEdgeOnSetOfGraphs,#,k]&//Select[#,(!GraphOfCutsHasBlownUpZeroPropQ[#])&]&//GenInternalLabelingsOfGraphOfCuts/@#&//Flatten;
);


(* ::Text:: *)
(*When writing up my code I should indicate who is memoized.  I want to memoize this function but if CubicBasisGraphs changes (because you look at a different amplitude without quitting the kernel) then this won't go well.  An alternative option is to make this a function of the CubicBasisGraphs.*)


CutLHS[GraphOfCuts_, AmpHead_]:=Module[{MomRules,dRules},
{MomRules,dRules}=GetGraphOfCutsKinematics[GraphOfCuts];
GraphOfCuts//Identity@@#&//DeleteCases[VerF[_]]//Times@@#&//#/.VerF[x__]:>AmpHead[{x}]&//#/.MomRules&//#/.dRules&
];


CutRHS[GraphOfCuts_,NumRepFunc_]:=Module[{MomCons,dRules,MomConsUncutMom,UncutMomenta,numer,denom},
Table[
{MomConsUncutMom,UncutMomenta}=GetBlownUpGraphMomConsAndUncutProps[GraphOfCuts,graph];
{MomCons,dRules}=GetGraphOfCutsKinematics[GraphOfCuts];
numer=graph//FindVerGraphNumerator;
denom=d2/@UncutMomenta//Times@@#&;
(numer//NumRepFunc)/denom/.MomConsUncutMom/.MomCons/.dRules
,{graph,GraphOfCuts//BlowUpOrderedCuts}]//Total
];


(* ::Subsection:: *)
(*NLSM 3-loop graphs*)


(* ::Input:: *)
(*(*graph1=Table[CubicBasisGraphs//Nest[PinchOneEdgeOnSetOfGraphs,#,k]&,{k,1,8(*6*)}]//Flatten;//AbsoluteTiming*)*)


(* ::Input:: *)
(*(*graph2=graph1//Select[#,(#//Identity@@#&//DeleteCases[VerF[_]]//#/.VerF[x__]:>Length[{x}]&//EvenQ/@#&//And@@#&)&]&//Select[#,FreeQ[#,VerF[___,x_,___,-x_,___]]&]&//Select[#,FreeQ[#,VerF[___,-x_,___,x_,___]]&]&;*)*)


(* ::Input:: *)
(*(*graph2//VerGraphHashCode/@#&*)*)


(* ::Subsection:: *)
(*End*)


End[];

EndPackage[];

Print["\n ---- PHYS GRAPH ---- \n\n The Phys Graph (or physics graph) package implements graph code for calculating scattering amplitudes.  The main purpose is to calculate BCJ representations of loop integrands so the pacakge includes code for graph symmetries (and isomorphisms), cuts, and Jacobi relations."];
