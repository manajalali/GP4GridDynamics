function Phi=Find_Diff3(t1,t2,beta,type,Kx)

switch type
case 00 
        Phi=Kx;
case 01
        Phi= 2*beta*(-Kx*diag(t2) + diag(t1)*Kx);
case 10 
        Phi= 2*beta*(Kx*diag(t2) - diag(t1)*Kx);
case 02
        Phi= -2*beta*Kx + (2*beta)^2*(diag(t1.*t1)*Kx +  Kx*diag(t2.*t2) - 2*diag(t1)*Kx*diag(t2));
case 20
        Phi= -2*beta*Kx + (2*beta)^2*(diag(t1.*t1)*Kx +  Kx*diag(t2.*t2) - 2*diag(t1)*Kx*diag(t2));
case 11
        Phi=  2*beta*Kx - (2*beta)^2*(diag(t1.*t1)*Kx +  Kx*diag(t2.*t2) - 2*diag(t1)*Kx*diag(t2));
case 12
        Phi= 3*(2*beta)^2*(diag(t1)*Kx - Kx*diag(t2));
        Phi= Phi - (2*beta)^3*(diag(t1.^3)*Kx - 3*diag(t1.^2)*Kx*diag(t2) + 3*diag(t1)*Kx*diag(t2.^2) - Kx*diag(t2.^3));
case 21
        Phi= -3*(2*beta)^2*(diag(t1)*Kx - Kx*diag(t2));
        Phi=Phi + (2*beta)^3*(diag(t1.^3)*Kx - 3*diag(t1.^2)*Kx*diag(t2) + 3*diag(t1)*Kx*diag(t2.^2) - Kx*diag(t2.^3));
case 22
        Phi=3*(2*beta)^2*Kx-6*(2*beta)^3*(diag(t1.*t1)*Kx +  Kx*diag(t2.*t2) - 2*diag(t1)*Kx*diag(t2));
        Phi=Phi+(2*beta)^4*(diag(t1.^4)*Kx-4*diag(t1.^3)*Kx*diag(t2)+6*diag(t1.^2)*Kx*diag(t2.^2)-4*diag(t1)*Kx*diag(t2.^3)+Kx*diag(t2.^4));     
end
end


