function [uL,uR] = F_VR1D_Reconstruct(xs,xc,rec,u,urec)

uL = u;
uR = u;

N = rec.N;

for icell = 1:N
    %left 
  
    DI_L = rec.DiBj_FL(1,:,icell);
    DI_R = rec.DiBj_FR(1,:,icell);
    
    
    UR = DI_R * urec(:,:,icell);
    UR(1,:) = UR(1,:)*1+u(:,icell)';
    UL = DI_L * urec(:,:,icell);
    UL(1,:) = UL(1,:)*1+u(:,icell)';
    
    uL(:,icell) = UL;
    uR(:,icell) = UR;
    
end