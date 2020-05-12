function result=phii(x)
    g=0.087;
    I=177.;
    c=615.; 
    y=c*x-I;
    result = y./(1-exp(-g*y));
end
