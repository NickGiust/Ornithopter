a = 3.5;
b = -.0334;
c = -.106;
d = .74;
m = 1.4;
n = 2;

Vb = @(s) a + b*(-log(s))^m + c*s + d*exp(n*(s-1));

fplot(Vb, [0 1]);