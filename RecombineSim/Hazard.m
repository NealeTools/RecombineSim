function y = Hazard(t,A,B)
A(A <= 0) = NaN;
B(B <= 0) = NaN;
y = gampdf(t,A,B)./(1-gamcdf(t,A,B));
y(t < 0) = 0;