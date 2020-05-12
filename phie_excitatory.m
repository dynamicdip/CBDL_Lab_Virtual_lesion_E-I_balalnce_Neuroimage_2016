function result=phie_excitatory(x)
a = 270;
b = 108;
d = 0.154;
y=a*x-b;
result = y./(1-exp(-d*y));
end
