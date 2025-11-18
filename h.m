function a = hopp(t)
    load(variabler.mat)
    if t <= L/v
        a = H/2*(1-cos(2*pi*v*t/L));
    else
        a = 0;
    end
end

hopp(0)

