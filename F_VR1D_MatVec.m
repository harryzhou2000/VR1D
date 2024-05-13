function [Murec, rhs] = F_VR1D_MatVec(urec,u,rec, WG)

N = rec.N;
nrec = rec.nrec;

rhs = urec * 0;
Murec = urec * 0;

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

    Murec(:,:,icell) = Aii * urec(:,:,icell);
    
    Murec(:,:,icell) = Murec(:,:,icell)-...
        ( Bij_L * urec(:,:,icellL) ) ;
    
    Murec(:,:,icell) = Murec(:,:,icell)-...
        ( Bij_R * urec(:,:,icellR) );

    rhs(:,:,icell) = rhs(:,:,icell) + ...
        ( b_L* ((u(:,icellL) - u(:,icell))') );
    rhs(:,:,icell) = rhs(:,:,icell) + ...
        ( b_R* ((u(:,icellR) - u(:,icell))') );
    
    
end

