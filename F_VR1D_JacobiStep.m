function urec = F_VR1D_JacobiStep(urec,u,rec, WG)

N = rec.N;
nrec = rec.nrec;

urecnew = urec * 0;


for i = 1:N
    icell = i;
    iFL = i;
    iFR = mod(i + 1 - 1, N) + 1;
    WGL = WG(iFL);
    WGR = WG(iFR);
    
    icellL = mod(icell - 1 - 1, N) + 1;
    icellR = mod(icell + 1 - 1, N) + 1;
    
    % % Face control
        Aii = rec.Aii_L(:,:,icell) * WGL + rec.Aii_R(:,:,icell) * WGR + rec.Aii_C(:,:,icell);
    % Vol control
%     Aii = rec.Aii_L(:,:,icell) * 1 + rec.Aii_R(:,:,icell) * 1 + rec.ALimit(:,:,icell) * max(0,-log(0.5*(WGL + WGR))) * 0;
    
    Bij_L =  rec.Bij_L(:,:,icell) * WGL;
    Bij_R =  rec.Bij_R(:,:,icell) * WGR;
    
    b_L =  rec.b_L(:,icell) * WGL;
    b_R =  rec.b_R(:,icell) * WGR;
    
    urecnew(:,:,icell) = urecnew(:,:,icell)+...
        Aii\( Bij_L * urec(:,:,icellL) ) + ...
        Aii\( b_L* ((u(:,icellL) - u(:,icell))') );
    
    urecnew(:,:,icell) = urecnew(:,:,icell)+...
        Aii\( Bij_R * urec(:,:,icellR) ) + ...
        Aii\( b_R* ((u(:,icellR) - u(:,icell))') );
    
end

urec = urecnew;