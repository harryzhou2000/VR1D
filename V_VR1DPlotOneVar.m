function [f] = V_VR1DPlotOneVar(f, xs, xc, u, urec,rec, ivar, nsamp, baseFuncType)

if(~exist("baseFuncType", 'var'))
    baseFuncType = 0;
end

N = size(xs,2)-1;
nrec = rec.nrec;

for icell = 1:N
    coords = xs(icell:icell+1);
    xC = xc(icell);
    xsamp = linspace(coords(1),coords(2),nsamp+1);
    xisamp = linspace(-1,1,nsamp+1);
    vsamp = 0 * xsamp;



    if baseFuncType == 0
        orthCoef = eye(rec.nrec);
        if(rec.useOrth)
            orthCoef = rec.orthCoef(:,:,icell);
        end
        for isamp = 1:numel(xsamp)
            vsamp(isamp) = ...
                F_VR1D_DiffBaseValue(xC,xisamp(isamp),coords,rec.baseMoment(:,icell)',zeros(1,nrec)) * orthCoef' * urec(1:nrec,ivar,icell) + ...
                u(ivar,icell);
        end
    else
        for isamp = 1:numel(xsamp)
            vsamp(isamp) = u(ivar,icell) + (rec.fphi(xsamp(isamp)) - rec.meanPhi(:,icell))' * urec(1:nrec,ivar,icell);
        end
    end
    hold on;
    plot(f,xsamp,vsamp);
    plot(xC,u(ivar,icell),'.');


end

hold off;