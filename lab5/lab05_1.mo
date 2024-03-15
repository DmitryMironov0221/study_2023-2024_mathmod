model lab05_1
Real a = 0.63;
Real b = 0.019;
Real c = 0.59;
Real d = 0.018;
Real x;
Real y;
initial equation
x = 7;
y = 12;
equation
der(x) = -a*x + b*x*y;
der(y) = c*y - d*x*y;
end lab05_1;