function result=phie(x)
g=0.16;
I=125.;
c=310.;
y=c*x-I;
result = y./(1-exp(-g*y));
end
