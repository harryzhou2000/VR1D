function v = F_VR1D_DiffBaseValue(xC, xparam, coords, moment, v)

reflen = (coords(2) - coords(1)) * 0.5;
x = xC + xparam * reflen;
xR = (x - xC)/reflen;


for idiff = 1:size(v,1)
    for ibase = 1:size(v,2)
        dx = idiff-1;
        px = ibase;
        c = 0;
        if(dx <= px)
            c = factorial(px)/factorial(px-dx);
        end
        v(idiff,ibase) = c * xR^(px-dx) / reflen^dx;
        if c==0
            v(idiff,ibase) = 0;
        end
    end
end

v(1,:) = v(1,:) - moment;