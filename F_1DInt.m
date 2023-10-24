function result = F_1DInt(f)
%  f(iint, xi)


xis = [-0.906179845938664, -0.538469310105683, 0, 0.538469310105683, 0.906179845938664];
wes = [0.236926885056190, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056190];

result = f(1,xis(1)) * wes(1);

for i = 2:numel(xis)
   result = result + f(i,xis(i)) * wes(i); 
end