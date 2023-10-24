function rec = F_VR1D_GetRecMat(xs,xc, ppoly,rec)

N = size(xs,2)-1;
nrec = ppoly;
WD = [1,1,1/2,1/6];
wGG = 1;

rec.nrec = nrec;
rec.N = N;

WD = WD(1:nrec+1);

rec.baseMoment = zeros(nrec,N);
rec.DiBj_FL = zeros(nrec+1,nrec,N);
rec.DiBj_FR = zeros(nrec+1,nrec,N);
rec.ALimit = zeros(nrec,nrec,N);
rec.Aii_C = zeros(nrec,nrec,N);
rec.Aii_L = zeros(nrec,nrec,N);
rec.Aii_R = zeros(nrec,nrec,N);
rec.Bij_L = zeros(nrec,nrec,N);
rec.Bij_R = zeros(nrec,nrec,N);
rec.b_L = zeros(nrec,N);
rec.b_R = zeros(nrec,N);

for icell = 1:N
    coords = xs(icell:icell+1);
    xC = xc(icell);
    baseMoment = F_1DInt(@(i,xi) F_VR1D_DiffBaseValue(xC,xi,coords,zeros(1,nrec),zeros(1,nrec)) );
    rec.baseMoment(:,icell) = baseMoment';
    
end

for icell = 1:N
    %left
    moment = rec.baseMoment(:,icell)';
    coords = xs(icell:icell+1);
    reflenI = coords(2) - coords(1);
    xC = xc(icell);
    
    icellL = mod(icell - 1 - 1, N) + 1;
    icellR = mod(icell + 1 - 1, N) + 1;
    coordsL = xs(icellL:icellL+1);
    coordsR = xs(icellR:icellR+1);
    xCL = xc(icellL);
    xCR = xc(icellR);
    momentL = rec.baseMoment(:,icellL)';
    momentR = rec.baseMoment(:,icellR)';
    
    reflenL = 0.5 * (coordsL(2)-coordsL(1) + coords(2) - coords(1));
    reflenR = 0.5 * (coordsR(2)-coordsR(1) + coords(2) - coords(1));
    WDGL = WD .* reflenL.^((1:nrec+1) -1);
    WDGR = WD .* reflenR.^((1:nrec+1) -1);
    WDGI = WD .* reflenI.^((1:nrec+1) -1);
    
    DI_L = F_VR1D_DiffBaseValue(xC,-1,coords,moment,zeros(nrec + 1, nrec));
    DI_R = F_VR1D_DiffBaseValue(xC, 1,coords,moment,zeros(nrec + 1, nrec));
    rec.DiBj_FL(:,:,icell) = DI_L;
    rec.DiBj_FR(:,:,icell) = DI_R;
    
    rec.Aii_L(:,:,icell) =  F_VR1D_FaceFunctional(DI_L,DI_L,WDGL);
    rec.Aii_R(:,:,icell) =  F_VR1D_FaceFunctional(DI_R,DI_R,WDGR);

    DL_R = F_VR1D_DiffBaseValue(xCL, 1,coordsL,momentL,zeros(nrec + 1, nrec));
    DR_L = F_VR1D_DiffBaseValue(xCR,-1,coordsR,momentR,zeros(nrec + 1, nrec));

    
    
    
   
    
    rec.Bij_L(:,:,icell) = F_VR1D_FaceFunctional(DI_L,DL_R,WDGL);
    rec.Bij_R(:,:,icell) = F_VR1D_FaceFunctional(DI_R,DR_L,WDGR);
    
    rec.b_L(:,icell) = F_VR1D_FaceFunctional(1,DI_L(1,:), WDGL(1))';
    rec.b_R(:,icell) = F_VR1D_FaceFunctional(1,DI_R(1,:), WDGR(1))';

    intDiBj = F_1DInt(@(iG, xi) F_VR1D_DiffBaseValue(xC,xi,coords,moment,zeros(nrec + 1, nrec)) );
    intD1Bj = intDiBj(2, :) * reflenI;
    A_GG_h = -(DI_R(1,:) - DI_L(1,:)) / 2 + intD1Bj;
    A_GG = A_GG_h' * A_GG_h;
    BR_GG_h = +DR_L(1,:) / 2;
    BL_GG_h = -DL_R(1,:) / 2;

    BR_GG = A_GG_h' * BR_GG_h;
    BL_GG = A_GG_h' * BL_GG_h;

    br_GG = A_GG_h' * ( 1/2);
    bl_GG = A_GG_h' * (-1/2);

    
    rec.Aii_C(:,:,icell) = rec.Aii_C(:,:,icell) +  A_GG * wGG;
    rec.Bij_L(:,:,icell) = rec.Bij_L(:,:,icell) + BL_GG * wGG;
    rec.Bij_R(:,:,icell) = rec.Bij_R(:,:,icell) + BR_GG * wGG;
    rec.b_L(:,icell) = rec.b_L(:,icell) + bl_GG * wGG;
    rec.b_R(:,icell) = rec.b_R(:,icell) + br_GG * wGG;

    ALimit = F_1DInt(@(i,xi) F_VR1D_FaceFunctional(...
        F_VR1D_DiffBaseValue(xC,xi,coords, moment,zeros(1 + nrec,nrec) ),...
        F_VR1D_DiffBaseValue(xC,xi,coords, moment,zeros(1 + nrec,nrec) ),...
        WDGI)...
        );
    rec.ALimit(:,:,icell) = ALimit;
    %     rec.Aii_L(:,:,icell)
end

%     function dbv = getDBVVol(xC,xi,coords,moment)
%         dbv = F_VR1D_DiffBaseValue(xC,xi,coords,moment,zeros(1 + nrec,nrec));
%         
%     end
% 
%     function A = ARegVol(xC,xi,coords,moment, WDGI)
%         dbv = getDBVVol(xC,xi,coords,moment);
%         A = F_VR1D_FaceFunctional(dbv,dbv, WDGI)';
%     end