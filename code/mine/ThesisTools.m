BeginPackage["ThesisTools`"]

GetKerrgll::usage="GetKerrgll[M,a,r,T] returns the Kerr metric in Boyer-Lindquist coordinates using the labels r and T for the radial and angular coordinate, M for mass, and a for spin."

RKODE4::usage="RKODE4[G,x0,f0,dx,Ns] is a 4th order Runge-Kutta integrator.  Refer to Kerr.Geodesics.nb for usage."

OuterErg::usage = "OuterErg[m,a,T] returns the radius of the outer ergosphere as a function of theta."

InnerErg::usage = "InnerErg[m,a,T] returns the radius of the inner ergosphere as a function of theta."

OuterEH::usage = "OuterErg[m,a,T] returns the radius of the outer event horizon as a function of theta."

InnerEH::usage = "InnerEH[m,a,T] returns the radius of the inner event horizon as a function of theta."

Surfaces::usage = "Surfaces[m,a,T] returns the radii of all the Kerr geometry surfaces in order of their appearance (Outer ergosphere, outer event horizon, inner event horizon, inner ergosphere)"

Begin["`Private`"]

GetKerrgll[M_, a_, r_, T_] := 
 FullSimplify[{{-1 + 2 M r/(a^2 Cos[T]^2 + r^2), 0, 0, -2 a M r Sin[T]^2/(r^2 + a^2 Cos[T]^2)}, {0, (a^2 Cos[T]^2 + r^2)/(a^2 - 2 M r + r^2), 0, 0}, {0, 0, a^2 Cos[T]^2 + r^2, 0}, {-2 a M r Sin[T]^2/(r^2 + a^2 Cos[T]^2), 0, 0, (Sin[T]^2 ((a^2 + r^2)^2 - a^2 (a^2 - 2 M r + r^2) Sin[T]^2))/(a^2 Cos[T]^2 + r^2)}}]

RKODE4[G_, x0_, f0_, dx_, Ns_] := 
 Module[{k1, k2, k3, k4, retvals, index, xn, fn},
  retvals = Table[{0.0, 0.0}, {i, 1, Ns}];
  xn = x0;
  fn = f0;
  retvals[[1]] = {x0, f0};
  For[index = 2, index <= Ns, index = index + 1,
   k1 = dx G[xn, fn];
   k2 = dx G[xn + dx/2, fn + k1/2];
   k3 = dx G[xn + dx/2, fn + k2/2];
   k4 = dx G[xn + dx, fn + k3];
   fn = fn + (1/3) (k1/2 + k2 + k3 + k4/2);
   xn = xn + dx;
   retvals[[index]] = {xn, fn};
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