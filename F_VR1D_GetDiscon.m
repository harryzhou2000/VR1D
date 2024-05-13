function discon = F_VR1D_GetDiscon(xs,xc,u,urec, rec)

N = rec.N;
discon = zeros(1,N);
WD = [1,1,1/2,1/6];
nrec = rec.nrec;
WD = WD(1:nrec+1);

for icell = 1:N
    %left
    coords = xs(icell:icell+1);
    
    icellL = mod(icell - 1 - 1, N) + 1;
    coordsL = xs(icellL:icellL+1);
    reflenL = 0.5 * (coordsL(2)-coordsL(1) + coords(2) - coords(1));
    
    WDGL = WD .* reflenL.^((1:nrec+1) -1);
    
    DI_L = rec.DiBj_FL(:,:,icell);
    DL_R = rec.DiBj_FR(:,:,icellL);
    
    
    UL = DL_R * urec(1:nrec,:,icellL);
    UL(1,:) = UL(1,:)+u(:,icellL)';
    UR = DI_L * urec(1:nrec,:,icell);
    UR(1,:) = UR(1,:)+u(:,icell)';
    
    UDIFF = UL-UR;
    USUM = UL+UR;
    IJI = F_VR1D_FaceFunctional(UDIFF,UDIFF,WDGL);
    ISI = F_VR1D_FaceFunctional(UDIFF,UDIFF,WDGL);
    
    
    discon(icell) = max(IJI,[],'all');
%     rec.Aii_L(:,:,icell)
end