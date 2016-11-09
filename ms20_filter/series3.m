function [b1, b2, b3] = series3(a1, Rp1, a2, Rp2, a3, Rp3)
% whatever...
    r = 2 * (a1 + a2 + a3) / (Rp1 + Rp2 + Rp3);
    
    b1 = a1 - Rp1*r;
	b2 = a2 - Rp2*r;
	b3 = a3 - Rp3*r;
end