function urec = F_VR1D_ArrayInit(u, rec)

nvar = size(u,1);
N = size(u,2);

nRecDOF = rec.ppoly;
if rec.wGGLM
    nRecDOF = nRecDOF + 1;
end
urec = zeros(nRecDOF,nvar,N);