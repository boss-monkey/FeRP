function lnphi = cal_lnphi(A,B,Z)
        sqrt2 = sqrt(2);
        E =  (Z + (1 + sqrt2).*B)./(Z + (1 - sqrt2).*B);
        lnphi = (Z - 1) - log(Z - B) - A./(2*sqrt2.*B).*log(E);
end
