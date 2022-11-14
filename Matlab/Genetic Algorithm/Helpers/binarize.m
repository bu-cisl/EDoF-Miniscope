function a = binarize(b,thres)
    temp = real(b);
    temp(temp >= thres) = max(real(b(:)));
    temp(temp < thres) = min(real(b(:)));
    a = temp;
end