function [urec,WG_dist] = F_VR1D_StaticRec_Simple(urec, u, rec, xs, xc, bcType)
if(~exist("bcType", 'var'))
    bcType = 0;
end
N = rec.N;
urec_lin = urec;

resTh_rec = 1e-8;
max_iter = 5000;

WG_lin = ones(1,N);
WG_dist = WG_lin;


%% use jacobi
for iter = 1:max_iter
    urec_lin = F_VR1D_JacobiStep(urec_lin, u, rec, WG_lin,bcType);
    urecNew = urec_lin;
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
%% use gmres
% [~,rhs] = F_VR1D_MatVec(urec_lin, u, rec, WG_lin);
% 
% [urecSol,~,~,iterGMRES, resGMRES ] = ...
%     gmres(@(urecc) reshape(F_VR1D_MatVec(reshape(urecc,size(urec)), u, rec, WG_lin),[],1),...
%     rhs(:), 20, resTh_rec, max_iter);
% urec = reshape(urecSol,size(urec));
% fprintf("resGMRES:\n");
% fprintf("%.4g\n", resGMRES);

