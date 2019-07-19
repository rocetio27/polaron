%%fourier transform
function val=fourier_transform(w,Nt,t,Phif,N)
    eV=1/27.211385;
    integrand = zeros(1,Nt);
    for i=1:Nt
            integrand(i) = -1i*exp(1i*w*t(i))*exp(-Phif(i))*exp(-1e-5*t(i)); %*exp(-1e-5*t(i)) *exp(-1i*0.01*eV*t(i))
    end
    dt=t(2)-t(1);
    val=0;
    for i=1:N
        val=val+2*dt/45*(7*integrand(5*(i-1)+1)+32*integrand(5*(i-1)+2)+12*integrand(5*(i-1)+3)+32*integrand(5*(i-1)+4)+7*integrand(5*(i-1)+5));
    end
end

% for i=1:6*N+3
%     for j=1:2*N+1
%         G(i)=G(i)+dt/3*(Integrand(i,1+3*(j-1))+4*Integrand(i,2+3*(j-1))+Integrand(i,3+3*(j-1)));
%     end
% end

% for i=1:6*N+3
%     for j=1:6*N
%         G(i)=G(i)+dt/2*(Integrand(i,j)+Integrand(i,j+1));
%     end
% end