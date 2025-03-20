function ans_s = sgn_BW(x)
    ans_s(x>=0) = 1;
    ans_s(x<0) = -1;
end