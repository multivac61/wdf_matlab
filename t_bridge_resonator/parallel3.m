function [b1, b2, b3] = parallel3(a1, Rp1, a2, Rp2, a3, Rp3)
% adapted for port 1
    Gp1 = 1/Rp1;
	Gp2 = 1/Rp2;
	Gp3 = 1/Rp3;
    R = 2*(Gp1*a1 + Gp2*a2 + Gp3*a3)/(Gp1 + Gp2 + Gp3);
    
    b1 = R - a1;
	b2 = R - a2;
	b3 = R - a3;
end