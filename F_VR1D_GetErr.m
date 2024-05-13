function err = F_VR1D_GetErr(xs, xc, u, urec, rec, ivar, fAna, baseFuncType)
if(~exist("baseFuncType", 'var'))
    baseFuncType = 0;
end

N = size(xs,2)-1;
nrec = rec.nrec;
err = 0;
for icell = 1:N
    coords = xs(icell:icell+1);
    h = (coords(2)  -coords(1));

    if(baseFuncType == 0)
        orthCoef = eye(rec.nrec);
        if(rec.useOrth)
            orthCoef = rec.orthCoef(:,:,icell);
        end
        errInt = F_1DInt(@(iG, xi)...
            abs(F_VR1D_DiffBaseValue(xc(icell),xi,coords,rec.baseMoment(:,icell)',zeros(1,nrec)) * orthCoef' * urec(1:nrec,ivar,icell) + ...
            u(ivar,icell) - fAna(coords(1) + h * (xi + 1)/2) ) ) / 2 * h;
    else
        errInt = F_1DInt(@(iG, xi)...
            abs((rec.fphi((xi+1)/2*h + coords(1)) - rec.meanPhi(:,icell))' * urec(1:nrec,ivar,icell) + ...
             u(ivar,icell) - fAna(coords(1) + h * (xi + 1)/2) ) ) / 2 * h;

    end


    err = err+ errInt;
end