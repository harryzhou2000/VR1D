function rec = F_VR1D_GetRecMat(xs,xc, ppoly,rec)

N = size(xs,2)-1;
nrec = ppoly;
WD = [1,1,1/2,1/6] * 0;
wGG = 0;
wGG_Beta = 0;
wGGLM = 0;
rec.useOrth=false;
% rec.useOrth = true;

%PVR
wProjBase = [1,1/1,1/5,1/32] * 1;
wProj = 1;

%PVR

rec.nrec = nrec;
rec.N = N;

WD = WD(1:nrec+1);

rec.baseMoment = zeros(nrec,N);
rec.DiBj_FL = zeros(nrec+1,nrec,N);
rec.DiBj_FR = zeros(nrec+1,nrec,N);
nExtra = 0;
if wGGLM
    nExtra = 1;
end
nrecE = nrec + nExtra;
rec.ALimit = zeros(nrecE,nrecE,N);
rec.Aii_C = zeros(nrecE,nrecE,N);
rec.Aii_L = zeros(nrecE,nrecE,N);
rec.Aii_R = zeros(nrecE,nrecE,N);
rec.Bij_L = zeros(nrecE,nrecE,N);
rec.Bij_R = zeros(nrecE,nrecE,N);
rec.b_L = zeros(nrecE,N);
rec.b_R = zeros(nrecE,N);
rec.ppoly = ppoly;
rec.wGGLM = wGGLM;

rec.orthCoef = zeros(nrec, nrec,N);
rec.gBase = zeros(nrec, nrec,N);

for icell = 1:N
    coords = xs(icell:icell+1);
    xC = xc(icell);
    baseMoment = F_1DInt(@(i,xi) F_VR1D_DiffBaseValue(xC,xi,coords,zeros(1,nrec),zeros(1,nrec)) ) / 2;
    rec.baseMoment(:,icell) = baseMoment';

    gBaseIJ = F_1DInt(@(i,xi) F_VR1D_DiffBaseValue(xC,xi,coords,baseMoment,zeros(1,nrec))' * F_VR1D_DiffBaseValue(xC,xi,coords,baseMoment,zeros(1,nrec)) ) / 2;
    orthCoef = inv(chol(gBaseIJ)');
    rec.orthCoef(:,:, icell) = (orthCoef);
    gBaseIJOrth = F_1DInt(@(i,xi) orthCoef * F_VR1D_DiffBaseValue(xC,xi,coords,baseMoment,zeros(1,nrec))' * F_VR1D_DiffBaseValue(xC,xi,coords,baseMoment,zeros(1,nrec)) * orthCoef' ) / 2;
end

for icell = 1:N
    orthCoef = eye(nrec);
    if rec.useOrth
        orthCoef = rec.orthCoef(:,:,icell);
    end
    gBase = F_1DInt(@(i,xi) orthCoef * F_VR1D_DiffBaseValue(xC,xi,coords,baseMoment,zeros(1,nrec))' * F_VR1D_DiffBaseValue(xC,xi,coords,baseMoment,zeros(1,nrec)) * orthCoef' ) / 2;
    rec.gBase(:,:,icell) = gBase;
end

periodicLen = max(xs,[],"all") - min(xs,[],"all");

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
    if(icell == 1)
        xCL = xCL - periodicLen;
        coordsL = coordsL - periodicLen;
    end
    if(icell == N)
        xCR = xCR + periodicLen;
        coordsR = coordsR + periodicLen;
    end
    JL = (coordsL(2) - coordsL(1)) / 2;
    JR = (coordsR(2) - coordsR(1)) / 2;
    J = (coords(2) - coords(1)) / 2;
    momentL = rec.baseMoment(:,icellL)';
    momentR = rec.baseMoment(:,icellR)';
    orthCoef = eye(nrec);
    orthCoefL = eye(nrec);
    orthCoefR = eye(nrec);
    if rec.useOrth
        orthCoef = rec.orthCoef(:,:,icell);
        orthCoefL = rec.orthCoef(:,:,icellL);
        orthCoefR = rec.orthCoef(:,:,icellR);
    end
    gBase = rec.gBase(:,:,icell);
    gBaseL = rec.gBase(:,:,icellL);
    gBaseR = rec.gBase(:,:,icellR);


    reflenL = 0.5 * (coordsL(2)-coordsL(1) + coords(2) - coords(1));
    reflenR = 0.5 * (coordsR(2)-coordsR(1) + coords(2) - coords(1));
    WDGL = WD .* reflenL.^((1:nrec+1) -1);
    WDGR = WD .* reflenR.^((1:nrec+1) -1);
    WDGI = WD .* reflenI.^((1:nrec+1) -1);
    DI_L = F_VR1D_DiffBaseValue(xC,-1,coords,moment,zeros(nrec + 1, nrec)) * orthCoef';
    DI_R = F_VR1D_DiffBaseValue(xC, 1,coords,moment,zeros(nrec + 1, nrec)) * orthCoef';

    rec.DiBj_FL(:,:,icell) = DI_L;
    rec.DiBj_FR(:,:,icell) = DI_R;



    DL_R = F_VR1D_DiffBaseValue(xCL, 1,coordsL,momentL,zeros(nrec + 1, nrec)) * orthCoefL';
    DR_L = F_VR1D_DiffBaseValue(xCR,-1,coordsR,momentR,zeros(nrec + 1, nrec)) * orthCoefR';






    Bij_L = F_VR1D_FaceFunctional(DI_L,DL_R,WDGL);
    Bij_R = F_VR1D_FaceFunctional(DI_R,DR_L,WDGR);

    b_L = F_VR1D_FaceFunctional(1,DI_L(1,:), WDGL(1))';
    b_R = F_VR1D_FaceFunctional(1,DI_R(1,:), WDGR(1))';

    intDiBj = F_1DInt(@(iG, xi) F_VR1D_DiffBaseValue(xC,xi,coords,moment,zeros(nrec + 1, nrec))*orthCoef' );
    intD1Bj = intDiBj(2, :)  * (reflenI / 2);
    A_GG_h = (DI_R(1,:) - DI_L(1,:)) / 2 - intD1Bj + ...
        wGG_Beta * ( -DI_R(2,:) * reflenR - DI_L(2,:) * reflenL );
    A_GG = A_GG_h' * A_GG_h;
    BR_GG_h = +DR_L(1,:) / 2 + wGG_Beta * DR_L(2,:) * reflenR;
    BIR_GG_h = -DI_R(1,:) / 2 + wGG_Beta * DI_R(2,:) * reflenR;
    BL_GG_h = -DL_R(1,:) / 2 + wGG_Beta * DL_R(2,:) * reflenL;
    BIL_GG_h = +DI_L(1,:) / 2 + wGG_Beta * DI_L(2,:) * reflenL;


    BL_GG = -A_GG_h' * BL_GG_h;
    BR_GG = -A_GG_h' * BR_GG_h;



    bl_GG = -A_GG_h' * (-1/2);
    br_GG = -A_GG_h' * ( 1/2);

    f_dbv1_L = @(xi) F_VR1D_DiffBaseValue(xC, (xi*JL + xCL -xC) / J,coords,moment,zeros(1, nrec))*orthCoef';
    f_dbv1_R = @(xi) F_VR1D_DiffBaseValue(xC, (xi*JR + xCR -xC) / J,coords,moment,zeros(1, nrec))*orthCoef';

    f_dbv1_LI = @(xi) F_VR1D_DiffBaseValue(xCL, (xi*J + xC - xCL) / JL,coordsL,momentL,zeros(1, nrec))* orthCoefL';
    f_dbv1_RI = @(xi) F_VR1D_DiffBaseValue(xCR, (xi*J + xC - xCR) / JR,coordsR,momentR,zeros(1, nrec))* orthCoefR';

    % PVR:
    bij_PVR_L = F_1DInt(@(iG, xi) f_dbv1_L(xi))'/2;
    bij_PVR_R = F_1DInt(@(iG, xi) f_dbv1_R(xi))'/2;

    BBij_PVR_IL = F_1DInt(@(iG, xi) f_dbv1_L(xi)' * F_VR1D_DiffBaseValue(xCL, xi,coordsL,momentL,zeros(1, nrec)) * orthCoefL' )/2;
    BBij_PVR_IR = F_1DInt(@(iG, xi) f_dbv1_R(xi)' * F_VR1D_DiffBaseValue(xCR, xi,coordsR,momentR,zeros(1, nrec)) * orthCoefR' )/2;

    BBij_PVR_LI = F_1DInt(@(iG, xi) f_dbv1_LI(xi)' * F_VR1D_DiffBaseValue(xC, xi,coords,moment,zeros(1, nrec)) * orthCoef')/2;
    BBij_PVR_RI = F_1DInt(@(iG, xi) f_dbv1_RI(xi)' * F_VR1D_DiffBaseValue(xC, xi,coords,moment,zeros(1, nrec)) * orthCoef')/2;

    WD_PVR = diag(wProjBase(2:1+nrec));
    A_PVR_L = gBase * WD_PVR * gBase + BBij_PVR_IL * WD_PVR * BBij_PVR_IL' + bij_PVR_L * bij_PVR_L' * wProjBase(1);
    A_PVR_R = gBase * WD_PVR * gBase + BBij_PVR_IR * WD_PVR * BBij_PVR_IR' + bij_PVR_R * bij_PVR_R' * wProjBase(1);
    B_PVR_L = BBij_PVR_IL * WD_PVR * gBaseL + gBase * WD_PVR * BBij_PVR_LI';
    B_PVR_R = BBij_PVR_IR * WD_PVR * gBaseR + gBase * WD_PVR * BBij_PVR_RI';

    Aii_C = A_GG * wGG;
    Bij_L = Bij_L + BL_GG * wGG + B_PVR_L;
    Bij_R = Bij_R + BR_GG * wGG + B_PVR_R;
    b_L = b_L + bl_GG * wGG + bij_PVR_L * wProjBase(1);
    b_R = b_R + br_GG * wGG + bij_PVR_R * wProjBase(1);

    if wGGLM
        Aii_C = [Aii_C, -A_GG_h' * wGGLM; -A_GG_h * wGGLM, 0];
        Bij_L = [Bij_L, BIL_GG_h'* wGGLM; BL_GG_h* wGGLM, 0];
        Bij_R = [Bij_R, BIR_GG_h'* wGGLM; BR_GG_h* wGGLM, 0];
        b_L = [b_L;-1/2* wGGLM];
        b_R = [b_R; 1/2* wGGLM];
    end

    rec.Aii_L(1:nrec,1:nrec,icell) =  F_VR1D_FaceFunctional(DI_L,DI_L,WDGL) + A_PVR_L;
    rec.Aii_R(1:nrec,1:nrec,icell) =  F_VR1D_FaceFunctional(DI_R,DI_R,WDGR) + A_PVR_R;
    rec.Aii_C(:,:,icell) = Aii_C;
    rec.Bij_L(:,:,icell) = Bij_L;
    rec.Bij_R(:,:,icell) = Bij_R;
    rec.b_L(:,icell) = b_L;
    rec.b_R(:,icell) = b_R;

    ALimit = F_1DInt(@(i,xi) F_VR1D_FaceFunctional(...
        F_VR1D_DiffBaseValue(xC,xi,coords, moment,zeros(1 + nrec,nrec) ) * orthCoef',...
        F_VR1D_DiffBaseValue(xC,xi,coords, moment,zeros(1 + nrec,nrec) ) * orthCoef',...
        WDGI)...
        );
    rec.ALimit(1:nrec,1:nrec,icell) = ALimit;
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