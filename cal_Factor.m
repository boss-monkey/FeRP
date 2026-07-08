function [ dF1,dF2,dF3 ] = cal_Factor(A,B,Z)

    dF1 = 3*Z.^2 -2*(1-B).*Z +(A-2*B-3*B.^2);
    dF2 = Z - B;
    dF3 = Z.^2 - 2*Z - 6*B.*Z -A +2.*B +3.*B.^2;     
    
end

