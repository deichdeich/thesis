(* ::Package:: *)

BeginPackage["ThesisTools`"]
DoThesis::usage = "DoThesis[\!\(\*
StyleBox[\"title\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"author\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"shittiness\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"advisor\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"]\",\nFontFamily->\"Menlo\"]\) researches, compiles, and writes a Reed College thesis of quality \!\(\*
StyleBox[\"shittiness\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\).  If an \!\(\*
StyleBox[\"advisor\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\)object is specified, DoThesis will pass \!\(\*
StyleBox[\"advisor\",\nFontFamily->\"Times New Roman\",\nFontSlant->\"Italic\"]\) very rough work at the last minute.  Returns LaTeX document and lasting feelings of unworthiness and low self esteem."

SingleParticleNewtonianForce::usage = "SingleParticleNewtonianForce[{\!\(\*
StyleBox[\"array\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"}\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"i\",\nFontSlant->\"Italic\"]\)] returns the force on the \!\(\*SuperscriptBox[
StyleBox[\"i\",\nFontSlant->\"Italic\"], \(th\)]\) particle due to all other particles with a Newtonian \!\(\*SuperscriptBox[
StyleBox[\"r\",\nFontSlant->\"Italic\"], \(-2\)]\) potential.  \!\(\*
StyleBox[\"array\",\nFontSlant->\"Italic\"]\) is the state array, and is of the form {{{\!\(\*SubscriptBox[
StyleBox[\"m\",\nFontSlant->\"Italic\"], \(1\)]\)},{\!\(\*SubscriptBox[
StyleBox[\"x\",\nFontSlant->\"Italic\"], \(1\)]\),\!\(\*SubscriptBox[
StyleBox[\"y\",\nFontSlant->\"Italic\"], \(1\)]\)},{\!\(\*SubscriptBox[
StyleBox[OverscriptBox[\"x\", \".\"],\nFontSlant->\"Italic\"], \(1\)]\),\!\(\*SubscriptBox[
StyleBox[OverscriptBox[\"y\", \".\"],\nFontSlant->\"Italic\"], \(1\)]\)}},{{\!\(\*SubscriptBox[
StyleBox[\"m\",\nFontSlant->\"Italic\"], \(2\)]\)},{\!\(\*SubscriptBox[
StyleBox[\"x\",\nFontSlant->\"Italic\"], \(2\)]\),\!\(\*SubscriptBox[
StyleBox[\"y\",\nFontSlant->\"Italic\"], \(2\)]\)},{\!\(\*SubscriptBox[
StyleBox[OverscriptBox[\"x\", \".\"],\nFontSlant->\"Italic\"], \(2\)]\),\!\(\*SubscriptBox[
StyleBox[OverscriptBox[\"y\", \".\"],\nFontSlant->\"Italic\"], \(2\)]\)}},\[Ellipsis]}."

xy2rad::usage = "xy2rad[{x,y}] converts Cartesian {x,y} co\[ODoubleDot]rdinates to {r,\[Theta]} with careful treatment of vectors with the form {0,y}.  By convention, the zero-vector is preserved between Cartesian and polar."

G2xy::usage = "Converts polar derivative vector to {{xd,yd},{xdd,ydd}}"

xy2radMotion::usage = "xy2radMotion[{{\!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"y\",\nFontSlant->\"Italic\"]\)},{\!\(\*
StyleBox[OverscriptBox[
StyleBox[\"x\",\nFontSlant->\"Italic\"], \".\"],\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[OverscriptBox[\"y\", \".\"],\nFontSlant->\"Italic\"]\)}},{\!\(\*
StyleBox[OverscriptBox[\"x\", \"..\"],\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[OverscriptBox[\"y\", \"..\"],\nFontSlant->\"Italic\"]\)}] takes motion in Cartesian co\[ODoubleDot]rdinates and returns motion in polar, {{\!\(\*
StyleBox[OverscriptBox[\"r\", \".\"],\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[OverscriptBox[\"\[Theta]\", \".\"],\nFontSlant->\"Italic\"]\)},\!\(\*
StyleBox[OverscriptBox[
RowBox[{\"{\", \"r\"}], \"..\"],\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[OverscriptBox[\"\[Theta]\", \"..\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"}\",\nFontSlant->\"Italic\"]\)}."

PolarVectorSum::usage = "PolarVectorSum[{{\!\(\*SubscriptBox[\(\[Theta]\), \(1\)]\),\!\(\*SubscriptBox[\(r\), \(1\)]\)},{\!\(\*SubscriptBox[\(\[Theta]\), \(2\)]\),\!\(\*SubscriptBox[\(r\), \(2\)]\)},\[Ellipsis]}] adds vectors in polar form."

GetKerrgll::usage = "GetKerrgll[M,a,r,T] returs the Kerr metric in Boyer-Lindquist coordinates using the labels r and T for the radial and angular coordinate, M for mass, and a for spin."

RKODE4::usage="RKODE4[G,x0,f0,dx,Ns] is a 4th order Runge-Kutta integrator.  Refer to Kerr.Geodesics.nb for usage."

OuterErg::usage = "OuterErg[m,a,T] returns the radius of the outer ergosphere as a function of theta."

InnerErg::usage = "InnerErg[m,a,T] returns the radius of the inner ergosphere as a function of theta."

OuterEH::usage = "OuterErg[m,a,T] returns the radius of the outer event horizon as a function of theta."

InnerEH::usage = "InnerEH[m,a,T] returns the radius of the inner event horizon as a function of theta."

Surfaces::usage = "Surfaces[m,a,T] returns the radii of all the Kerr geometry surfaces in order of their appearance (Outer ergosphere, outer event horizon, inner event horizon, inner ergosphere)"

Begin["`Private`"]




SingleParticleNewtonianForce[f_,masses_,i_]:=Module[{j,x1,x2,y1,y2,m1,m2,distance2,jforce,jforcedir,totalforce,forces,Nparticles},
	Nparticles = Length[f];
	forces = Table[{0.,0.},Nparticles];
	x1 = f[[i,1,1]];
	y1 = f[[i,1,2]];
	m1 = masses[[i,1]];
	For [j=1,j<= Nparticles,j++,
		If[j!= i,
			m2 = masses[[j,1]];
			x2 = f[[j,1,1]];
			y2 = f[[j,1,2]];
			distance2 = ((x2-x1)^2+(y2-y1)^2);
			jforce = (m2)/distance2;
			jforcedir = {x2-x1,y2-y1}/Sqrt[distance2];
			forces[[j]] = jforce*jforcedir;
		];
	];
	Return[Total[forces]];
]
G2xy[G_,r_,T_]:=Module[{rd,rdd,Td,Tdd,xd,xdd,yd,ydd},
	rd = G[[1,1]];
	rdd = G[[1,2]];
	Td = G[[2,1]];
	Tdd = G[[2,2]];
	xd = Cos[T] rd-r Sin[T] Td;
	xdd = Cos[T] (-r Td^2+rdd)-Sin[T] (2 rd Td+r Tdd);
	yd = Sin[T] rd+Cos[T] r Td;
	ydd = 2 Cos[T] rd Td+Sin[T] (-r Td^2+rdd)+Cos[T] r Tdd;
	Return[{{xd,yd},{xdd,ydd}}];
]
xy2radMotion[state_,forces_]:=Module[{x,y,xd,yd,xdd,ydd,r2,r,rd,Td,rdd,Tdd},
	x = state[[1,1]];
	y = state[[1,2]];
	xd =state[[2,1]];
	yd = state[[2,2]];
	xdd = forces[[1]];
	ydd = forces[[2]];
	r2 = x^2+y^2;
	r = Sqrt[r2];
	rd = (x xd+y yd)/r;
	Td = (x yd-y xd)/r2;
	rdd = (-4 (x xd+y yd)^2+4 r2 (xd^2+yd^2+x xdd+y ydd))/(4 (r2)^(3/2));
	Tdd = (y^2 (-xdd y+2 xd yd)-x^2 (xdd y+2 xd yd)+x^3 ydd+x y (2 xd^2-2 yd^2+y ydd))/r2^2;
	Return[{{rd,Td},{rdd-r*Td^2,r*Tdd+2rd*Td}}];
]
xy2rad[point_] := Module[{r,phi,x,y},
	x = point[[1]];
	y = point[[2]];
	r = Sqrt[x^2+y^2];
	If[x>0,
		phi = ArcTan[y/x];
	];
	If[x<0,
		If[y>= 0,
			phi = ArcTan[y/x]+Pi;
		];
		If[y<0,
			phi = ArcTan[y/x]-Pi;
		];
	];
	If[x==0,
		If[y>0,
			phi = Pi/2;
		];
		If[y<0,
			phi = -Pi/2;
		];
		If[y==0,
			phi = 0;
		];
	];
	Return[{r,phi}];
]

PolarVectorSum[vectors_]:= Module[{v1,v2,r1,p1,r2,p2,r3,p3,i,Nvecs},
	Nvecs = Length[vectors];
	v1 = vectors[[1]];
	p1 =v1[[1]];
	r1 = v1[[2]];
	For[i = 2,i<=Nvecs,i++,
		v2 = vectors[[i]];
		r2 = v2[[2]];
		p2 = v2[[1]];
		r3 = Sqrt[r1^2 +r2^2+2 r1 r2 Cos[(p2-p1)]];
		p3 = p1+ArcCos[(r1+r2 Cos[(p2-p1)])/r3];
		r1 = r3;
		p1 = p3;
	];
	Return[{p1,r1}];
]

GetKerrgll[M_, a_, r_, T_] := 
 FullSimplify[{{-1 + 2 M r/(a^2 Cos[T]^2 + r^2), 0, 0, -2 a M r Sin[T]^2/(r^2 + a^2 Cos[T]^2)}, {0, (a^2 Cos[T]^2 + r^2)/(a^2 - 2 M r + r^2), 0, 0}, {0, 0, a^2 Cos[T]^2 + r^2, 0}, {-2 a M r Sin[T]^2/(r^2 + a^2 Cos[T]^2), 0, 0, (Sin[T]^2 ((a^2 + r^2)^2 - a^2 (a^2 - 2 M r + r^2) Sin[T]^2))/(a^2 Cos[T]^2 + r^2)}}]

RKODE4[G_, t0_, f0_, dt_, timesteps_] := 
 Module[{k1, k2, k3, k4, retvals, index, tn, fn},
  retvals = Table[{0.0, 0.0}, {i, 1, timesteps}];
  tn = t0;
  fn = f0;
  retvals[[1]] = {t0, f0};
  For[index = 2, index <= timesteps, index = index + 1,
   k1 = dt G[tn, fn];
   k2 = dt G[tn + dt/2, fn + k1/2];
   k3 = dt G[tn + dt/2, fn + k2/2];
   k4 = dt G[tn + dt, fn + k3];
   fn = fn + (1/3) (k1/2 + k2 + k3 + k4/2);
   tn = tn + dt;
   retvals[[index]] = {tn, fn};
   ];
  Return[retvals];
  ]

OuterErg [m_, a_, T_] := m + Sqrt[m^2 - a^2 Cos[T]^2];
InnerErg [m_, a_, T_] := m - Sqrt[m^2 - a^2 Cos[T]^2];
OuterEH[m_, a_, T_] := m + Sqrt[m^2 - a^2];
InnerEH[m_, a_, T_] := m - Sqrt[m^2 - a^2];

Surfaces[m_, a_, T_]:=FullSimplify[{OuterErg[m, a, T], OuterEH[m, a, T], InnerEH[m, a, T], InnerErg[m, a, T]}]

End[]
EndPackage[]






