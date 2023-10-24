function err = F_VR1D_GetErr(xs, xc, u, urec, rec, ivar, fAna)


N = size(xs,2)-1;
nrec = size(urec,1);
err = 0;
for icell = 1:N
    coords = xs(icell:icell+1);
    h = (coords(2)  -coords(1));
    
    errInt = F_1DInt(@(iG, xi)...
           abs(F_VR1D_DiffBaseValue(xc(icell),xi,coords,rec.baseMoment(:,icell)',zeros(1,nrec)) * urec(:,ivar,icell) + ...
           u(ivar,icell) - fAna(coords(1) + h * (xi + 1)/2) ) ) / 2 * h;
    
    
    
    err = err+ errInt;
end