function [f] = V_VR1DPlotOneVar(f, xs, xc, u, urec,rec, ivar, nsamp)

N = size(xs,2)-1;
nrec = size(urec,1);

for icell = 1:N
    coords = xs(icell:icell+1);
    xC = xc(icell);
    xsamp = linspace(coords(1),coords(2),nsamp+1);
    xisamp = linspace(-1,1,nsamp+1);
    vsamp = 0 * xsamp;
    
    for isamp = 1:numel(xsamp)
       vsamp(isamp) = ...
           F_VR1D_DiffBaseValue(xC,xisamp(isamp),coords,rec.baseMoment(:,icell)',zeros(1,nrec)) * urec(:,ivar,icell) + ...
           u(ivar,icell);
    end
    hold on;
    plot(f,xsamp,vsamp);
    plot(xC,u(ivar,icell),'.');
    
   
end

hold off;