clear;
nvar = 1;
ppoly = 3;

N = 100;
xs = linspace(-1,1,N+1) + 0;
xc = 0.5*(xs(2:end) + xs(1:end-1));

% u = sin(xc*2*pi);
fAna = @(xc) double(abs(xc - 0.5)<0.25) * (1) + cos(xc * 1* pi);
% fAna = @(xc) xc.^4;

u = fAna(xc);
for icell = 1:N
   u(:,icell) = F_1DInt(@(iG, xi) fAna(xs(icell) + (xs(icell+1)-xs(icell)) * (xi+1)/2)) / 2;
end

rec.name = "rec";
rec = F_VR1D_GetRecMat(xs,xc,ppoly,rec);
urec = F_VR1D_ArrayInit(u,rec);

recGCVR.name = "recGCVR";
recGCVR = F_GCVR1D_Init(xs,xc,ppoly,recGCVR);

%%
% clf;
% for iiter = 1:10
%     urecnew = F_VR1D_JacobiStep(urec,u,rec,ones(size(xc)));
%     inc = urecnew - urec;
%     urec = urecnew;
%     res = sum(abs(inc),'all');
%     fprintf('iter %d res %g\n', iiter, res);
%     cla;
%     
% end
% V_VR1DPlotOneVar(gca,xs,xc,u,urec,rec,1,10);
%     drawnow;
% discon = F_VR1D_GetDiscon(xs,xc,u,urec,rec);
% % discon = 1./(1+exp((discon-0.2)*32));
% discon = exp(-discon*32);
%%
% for iiter = 1:100
%     urecnew = F_VR1D_JacobiStep(urec,u,rec,discon);
%     inc = urecnew - urec;
%     urec = urecnew;
%     res = sum(abs(inc),'all');
%     fprintf('iter %d res %g\n', iiter, res);
%     cla;
% end
% V_VR1DPlotOneVar(gca,xs,xc,u,urec,rec,1,10);
% drawnow;

%% VR
cla;
% [urec, WG] = F_VR1D_StaticRec_C0(urec,u,rec,xs,xc);
[urec, WG] = F_VR1D_StaticRec_Simple(urec,u,rec,xs,xc,0);
V_VR1DPlotOneVar(gca,xs,xc,u,urec,rec,1,10,0);
err = F_VR1D_GetErr(xs, xc, u, urec, rec, 1, fAna);

fprintf("abs err = %e\n", err);

%% GCVR
% urec = urec * 0;
% cla;
% % [urec, WG] = F_VR1D_StaticRec_C0(urec,u,rec,xs,xc);
% [urec, WG] = F_VR1D_StaticRec_Simple(urec,u,recGCVR,xs,xc,0);
% V_VR1DPlotOneVar(gca,xs,xc,u,urec,recGCVR,1,10,1);
% err = F_VR1D_GetErr(xs, xc, u, urec, recGCVR, 1, fAna,1);
% 
% fprintf("abs err = %e\n", err);




