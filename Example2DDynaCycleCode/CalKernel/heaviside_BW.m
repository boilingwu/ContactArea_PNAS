function ans_h = heaviside_BW(x)
%     if x>0
%         ans = 1;
%     else 
%         ans = 0;
%     end
    ans_h(x>0) = 1;
    ans_h(x<=0) = 0;
end