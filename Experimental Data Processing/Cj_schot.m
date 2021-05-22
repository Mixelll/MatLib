%Plot capacitance of a schottky junction with doping N and area A
range = [-1 0.4];
N=5e16;
bl = 0.5;
R = 100e-4;
A = pi*R^2;
q = 1.6e-19;
es = 11.68 * 8.854e-14;
qes = q*es/2;
v = linspace(range(1), range(2));
c = A*sqrt(qes*N./abs(bl-v));
plot(v,c)

