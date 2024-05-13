function rec = F_GCVR1D_Init(xs, xc, ppoly,rec)


wd = [1e-1,1e-2,1e-3];
N = size(xs,2)-1;
fphi = @(x) (x.^(1:ppoly))';

nrec = ppoly;
wd = wd(1:nrec);


rec.nrec = nrec;
rec.N = N;
rec.fphi = fphi;

rec.meanPhi = nan(nrec, N);
rec.vol = xs(2:end) - xs(1:end-1);
rec.Aii_C = zeros(nrec, nrec, N);
rec.Aii_L = nan(nrec, nrec, N);
rec.Aii_R = nan(nrec, nrec, N);
rec.Bij_L = nan(nrec, nrec, N);
rec.Bij_R = nan(nrec, nrec, N);
rec.b_L = zeros(nrec, N);
rec.b_R = zeros(nrec, N);



for iC = 1:N
    x0 = xs(iC);
    x1 = xs(iC + 1);
    intBase = F_1DInt(@(i,xi) fphi(x0 + (xi + 1)/2 * (x1-x0))) / 2;
    rec.meanPhi(:, iC) = intBase;
    
end

for icell = 1:N
    icellL = mod(icell - 1 - 1, N) + 1;
    icellR = mod(icell + 1 - 1, N) + 1;

    volL = rec.vol(icellL);
    volR = rec.vol(icellR);
    volC = rec.vol(icell);

    AiiL = (volL + volC)^2 * ((rec.meanPhi(:, icell) * rec.meanPhi(:, icell)') + diag(wd));
    AiiR = (volR + volC)^2 * ((rec.meanPhi(:, icell) * rec.meanPhi(:, icell)') + diag(wd));
    BijL = (volL + volC)^2 * ((rec.meanPhi(:, icell) * rec.meanPhi(:, icellL)')+ diag(wd));
    BijR = (volR + volC)^2 * ((rec.meanPhi(:, icell) * rec.meanPhi(:, icellR)')+ diag(wd));
    bL = -(volL + volC)^2 * rec.meanPhi(:, icell);
    bR = -(volR + volC)^2 * rec.meanPhi(:, icell);

    rec.Aii_L(:,:,icell) = AiiL;
    rec.Aii_R(:,:,icell) = AiiR;
    rec.Bij_L(:,:,icell) = BijL;
    rec.Bij_R(:,:,icell) = BijR;
    rec.b_L(:,icell) = bL;
    rec.b_R(:,icell) = bR;

end



