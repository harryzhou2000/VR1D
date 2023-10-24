function [urec,WG_dist] = F_VR1D_StaticRec(urec, u, rec, xs, xc)

N = rec.N;
urec_lin = urec;

resTh_rec = 1e-4;
max_iter = 1000;

WG_lin = ones(1,N);


for iter = 1:max_iter
    urec_lin = F_VR1D_JacobiStep(urec_lin, u, rec, WG_lin);
    WG_dist = F_VR1D_GetDiscon(xs, xc, u, urec_lin, rec);
    WG_dist_L = circshift(WG_dist,1,2);
    WG_dist_R = circshift(WG_dist,-1,2);
%     WG_dist = max([WG_dist_L;WG_dist_R;WG_dist],[],1);
    WG_dist = exp(-WG_dist*8);
%     WG_dist(WG_dist<0.6) = WG_dist(WG_dist<0.6) * 1e-3;
    urecNew = F_VR1D_JacobiStep(urec, u, rec, WG_dist);
    inc = sum(abs(urecNew - urec),'all');
    urec = urecNew;
    
    if(iter == 1)
       inc0 = inc; 
    end
    fprintf("F_VR1D_StaticRec %d res = %g\n", iter, inc/inc0);
    if inc/inc0 < resTh_rec
       break; 
    end
end

