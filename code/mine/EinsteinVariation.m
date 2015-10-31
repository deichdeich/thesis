BeginPackage["EinsteinVariation`"]

MetricFromDS::usage=
"MetricFromDS[ds,X] takes a line element and generates
	a matrix form for g_mn given the coordinates X."

LineltFromMetric::usage="LineltFromMetric[gll,X] gives the line element associated with a metric g_mn and coordinates X - the returned form has simplified and grouped coefficients
as one would expect."

GetCCull::usage="GetCCull[gll,X] returns the Christoffel connection with lower symmetric indices and an upper index (Christoffel symbol of the second kind) constructed from the metric g_mn and the coordinates X."

GetCClll::usage="GetCClll[gll,X] returns the Christoffel connection for the metric g_mn and coordinates X, this is the Christoffel symbol of the first kind."

GetRiemann::usage="GetRiemann[gll,X] gives the Riemann tensor associated with the metric g_mn with coordinates X.  Note that the natural form is returned, with one upper, three lower indices: R^u_{ abg}."
GetKerrgll::usage="GetKerrgll[M,a,r,T] returns the Kerr metric in Boyer-Lindquist coordinates using the labels r and T for the radial and angular coordinate, M for mass, and a for spin."

GetRicci::usage="GetRicci[gll,X] computes the Ricci tensor (covariant form) for the metric g_mn with coordinates X."

GetRicciS::usage="GetRicciS[gll,X] computes the Ricci scalar for the metric g_mn with coordinates X."

Begin["`Private`"]

MetricFromDS[ds_,X_,simpopt_:{TimeConstraint->10}]:=Module[{dX,gll,iindex,jindex,dsi,Dim},Dim=Length[X];
gll=Table[0,{iXX,1,Dim},{jXX,1,Dim}];
dX=Table[ToExpression["d"<>ToString[X[[jXX]]]],{jXX,1,Dim}];
dsi=Expand[ds];
For[iindex=1,iindex<=Dim,iindex=iindex+1,gll[[iindex,iindex]]=Simplify[Coefficient[dsi,dX[[iindex]],2]];
For[jindex=iindex+1,jindex<=Dim,jindex=jindex+1,gll[[iindex,jindex]]=(1/2) Simplify[Coefficient[dsi,dX[[iindex]] dX[[jindex]],1]];
gll[[jindex,iindex]]=gll[[iindex,jindex]];];];
Return[Simplify[gll,simpopt]];]

LineltFromMetric[gll_,X_]:=Module[{dX,ds,Dim,iindex,jindex},Dim=Length[X];
dX=Table[ToExpression["d"<>ToString[X[[jXX]]]],{jXX,1,Dim}];
ds=Sum[Sum[If[jXX==kXX,Simplify[gll[[jXX,kXX]]] dX[[jXX]] dX[[kXX]],2 Simplify[gll[[jXX,kXX]]] dX[[jXX]] dX[[kXX]]],{kXX,jXX,Dim}],{jXX,1,Dim}];
Return[ds];]

GetCClll[gll_,X_,simpopt_:{TimeConstraint->10}]:=Module[{CClllout,Dim},Dim=Length[X];
CClllout=(1/2) Table[D[gll[[bXX,uXX]],X[[vXX]]]+D[gll[[bXX,vXX]],X[[uXX]]]-D[gll[[uXX,vXX]],X[[bXX]]],{bXX,1,Dim},{uXX,1,Dim},{vXX,1,Dim}];
Return[Simplify[CClllout,simpopt]];]

GetCCull[gll_,X_,simpopt_:{TimeConstraint->10}]:=Module[{CCullout,guu,CClllin,Dim},Dim=Length[X];
guu=Simplify[Inverse[gll]];
CClllin=GetCClll[gll,X];
CCullout=Table[Sum[guu[[aXX,bXX]] CClllin[[bXX,uXX,vXX]],{bXX,1,Dim}],{aXX,1,Dim},{uXX,1,Dim},{vXX,1,Dim}];
Return[Simplify[CCullout,simpopt]];]

GetRiemann[gll_,X_,simpopt_:{TimeConstraint->10}]:=Module[{Rulllout,CCull,Dim},Dim=Length[X];
CCull=GetCCull[gll,X,simpopt];
Rulllout=Table[D[CCull[[aXX,bXX,pXX]],X[[gXX]]]-D[CCull[[aXX,gXX,pXX]],X[[bXX]]]+Sum[CCull[[aXX,gXX,sXX]] CCull[[sXX,bXX,pXX]],{sXX,1,Dim}]-Sum[CCull[[aXX,bXX,sXX]] CCull[[sXX,gXX,pXX]],{sXX,1,Dim}],{aXX,1,Dim},{pXX,1,Dim},{gXX,1,Dim},{bXX,1,Dim}];
Return[Simplify[Rulllout,simpopt]];]

GetRicci[gll_,X_,simpopt_:{TimeConstraint->10}]:=Module[{Rulll,Rll,Dim},Dim=Length[X];
Rulll=GetRiemann[gll,X,simpopt];
Rll=Table[Sum[Rulll[[s1XX,aXX,s1XX,bXX]],{s1XX,1,Dim}],{aXX,1,Dim},{bXX,1,Dim}];
Return[Simplify[Rll,simpopt]];]

GetRicciS[gll_,X_,simpopt_:{TimeConstraint->10}]:=Module[{Rs,guu,Dim},Dim=Length[X];
Rll=GetRicci[gll,X,simpopt];
guu=Simplify[Inverse[gll],simpopt];
Rs=Sum[Rll[[s1XX,s2XX]] guu[[s1XX,s2XX]],{s1XX,1,Dim},{s2XX,1,Dim}];
Return[Simplify[Rs,simpopt]];]

GetKerrgll[M_,a_,r_,T_]:=FullSimplify[{{-1+2 M r/(a^2 Cos[T]^2+r^2),0,0,-2 a M r Sin[T]^2/(r^2+a^2 Cos[T]^2)},{0,(a^2 Cos[T]^2+r^2)/(a^2-2 M r+r^2),0,0},{0,0,a^2 Cos[T]^2+r^2,0},{-2 a M r Sin[T]^2/(r^2+a^2 Cos[T]^2),0,0,(Sin[T]^2 ((a^2+r^2)^2-a^2 (a^2-2 M r+r^2) Sin[T]^2))/(a^2 Cos[T]^2+r^2)}}]


End[]
EndPackage[]
