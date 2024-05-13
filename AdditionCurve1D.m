AC =[208.0000         0  224.0000
         0  152.8889         0
  224.0000         0  688.0000];

A_GG =[1.0000         0    1.0000
         0         0         0
    1.0000         0    1.0000];


ws = logspace(-5,5,10001);
cs = nan(size(ws));


for i = 1:numel(ws)
    cs(i) = cond(AC + ws(i) * A_GG);

end

%%
plot(ws,cs);
set(gca,"XScale", "log")
set(gca,"YScale", "log")

  