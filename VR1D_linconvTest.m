
clear;
CFL = 0.5;
tmax = 1;
see = 1;

nvar = 1;
ppoly = 3;

N = 128;
xs = linspace(0,1,N+1);
xc = 0.5*(xs(2:end) + xs(1:end-1));
a = 1;
h = xs(2) - xs(1);
hs = xs(2:end) - xs(1:end-1);
dt = CFL * h / a;
niter = round(tmax/dt);


% u = sin(xc*2*pi);
u = double(abs(xc - 0.5)<0.25) * (1) + cos(xc * 2* pi)*(-0);


field.xc = xc;
field.xs = xs;
field.hs = hs;
field.urec = F_VR1D_ArrayInit(u,ppoly);
field.rec = F_VR1D_GetRecMat(xs,xc,ppoly,[]);
field.ithis = 1:N;
field.ile = circshift(field.ithis, 1,2);
field.iri = circshift(field.ithis,-1,2);
field.a = a;

%%
u0 = u;
t = 0;
for iiter = 1:niter
    [dudt0,field] = frhs(u,field);
    u1 = u + dt * dudt0;
    cla;
    V_VR1DPlotOneVar(gca,field.xs,field.xc,u,field.urec,field.rec,1,10);
    xlim([0,1]);
    ylim('auto');
    
    [dudt1,field] = frhs(u1,field);
    unew = u + dt * 0.5 * dudt0 + dt * 0.5 * dudt1;
    
    u = unew;
    t = t + dt;
    if(mod(iiter,see) == 0 || iiter == niter)
        cla;
        V_VR1DPlotOneVar(gca,field.xs,field.xc,u,field.urec,field.rec,1,10);
        xlim([0,1]);
        ylim('auto');
        drawnow;
    end
    
end

%%
function [dudt,field] = frhs(u,field)

eps_lim = 10; %要足够大 1+
N_lim = 4; %~10

a = field.a;
[field.urec,field.WG] = F_VR1D_StaticRec_C0(field.urec, u, field.rec,field.xs,field.xc);
%     N = field.rec.N;
%     urecsiz = size(field.urec);
%     ppoly = field.rec.nrec;
%     field.urecLim = field.urec;

[field.uL,field.uR] = F_VR1D_Reconstruct(field.xs,field.xc,field.rec,u,field.urec);



fL_uL = field.uR(field.ile);
fL_uR = field.uL;

fL_flux = 0.5*a*(fL_uL + fL_uR) - 0.5*abs(a)*(fL_uR - fL_uL);

dudt = (fL_flux - fL_flux(field.iri)) ./ field.hs;

end


