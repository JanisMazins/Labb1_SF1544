load(variabler.mat)
function a = h(t)
    if t <= L/v
        a = H/2*(1-cos(2*pi*v*t/L));
    else
        a = 0;
    end
end