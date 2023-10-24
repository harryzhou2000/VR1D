
nvar = 1;
ppoly = 3;

N = 32;
xs = linspace(0,1,N+1);
xc = 0.5*(xs(2:end) + xs(1:end-1));

% u = sin(xc*2*pi);
u = double(abs(xc - 0.5)<0.25) * (0) + cos(xc * 2* pi);


urec = F_VR1D_ArrayInit(u,ppoly);

rec = [];
rec = F_VR1D_GetRecMat(xs,xc,ppoly,rec);
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

%%
cla;
[urec, WG] = F_VR1D_StaticRec_C0(urec,u,rec,xs,xc);
V_VR1DPlotOneVar(gca,xs,xc,u,urec,rec,1,10);








