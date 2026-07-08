function dlnphidT = cal_dlnphidT(A,B,Z,a,dadT,T)
    sqrt2 = 2^0.5;
    dAdT = A./a.*dadT - 2*A./T;
    dBdT = -B./T;
    dABdT = dAdT./B - dBdT.*A./B.^2;
    [ dFdZ,dFdA,dFdB ] = cal_Factor(A,B,Z);

    dZdT = -(dFdA.*dAdT+dFdB.*dBdT)./dFdZ;

    dlnEdT = 2*sqrt2*(Z.*dBdT - B.*(dZdT))./...
          ((Z+(1+sqrt2).*B).*(Z+(1-sqrt2).*B)); 
    E =  (Z + (1 + sqrt2).*B)./(Z + (1 - sqrt2).*B);
    dlnphidT = dZdT - 1./(Z-B).*(dZdT -dBdT ) - 1/(2*sqrt2)*(dABdT.*log(E) + A./B.*dlnEdT);


end

