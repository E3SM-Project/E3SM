function P = psteady(p,Po,vb,W)
% PSTEADY  Computes P(W) in steady state.  Vectorized in all arguments
% except p, the parameters structure.
% CODE WRITTEN BY ED BUELER: https://github.com/bueler/hydrolakes/tree/master/codes


if any(any(Po < 0))
  error('psteady() requires nonnegative overburden pressure Po'), end
if any(any(vb < 0))
  error('psteady() requires nonnegative sliding speed vb'), end

sbcube = p.c1 * vb / (p.c2 * p.A);
frac = max(0.0, p.Wr - W) ./ (W + p.Y0);
P = Po - (sbcube .* frac).^(1/3);
P(P < 0) = 0.0;

