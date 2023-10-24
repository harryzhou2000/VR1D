function urec = F_VR1D_ArrayInit(u, ppoly)

nvar = size(u,1);
N = size(u,2);

urec = zeros(ppoly,nvar,N);