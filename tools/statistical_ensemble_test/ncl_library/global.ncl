; Global variable, any stddev < tol is treated as 0.
; Setting it here ensures the same value is used for both creating ensembles
; and comparing a single run to the ensemble. 
tol = 1d-12

undef("sum_sqr")
function sum_sqr(x)
begin
  return(sum(x*x))
end
