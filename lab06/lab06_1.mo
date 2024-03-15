model lab06_1
Real N = 11700;
Real I;
Real R;
Real S;
Real alpha = 0.6;
Real beta = 0.2;
initial equation
I = 270;
R = 49;
S = N - I - R;
equation
der(S) = 0;
der(I) = -beta*I;
der(R) = beta*I;
end lab06_1;