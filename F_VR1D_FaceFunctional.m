function coj = F_VR1D_FaceFunctional(diffi,diffj,weights)
coj = diffi' * diag(weights.^2) * diffj;