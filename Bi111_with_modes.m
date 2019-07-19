tic


%% units
h=2*pi;
v=8e+5/2.1876912633e+6; %4.1e+5/2.1876912633e+6 a.u. SI/a.u. 
eV=1/27.211385; %a.u.
nm=1e-9/5.2917721092e-11; %a.u. SI/a.u.;
kb=8.617e-5*27.21; % in a.u.
T=0;
e=1; % SI/a.u.
k=0;
%% parameters
a=0.4538*nm;
Z=83;
M=3.4702009e-25/9.10938291e-31; %SI/a.u.;

%% Mq Define
mode_num=6;
Mq=zeros(mode_num,Na);
Mq_square=zeros(sqrt(Na),sqrt(Na),mode_num);
for mode=1:6
    for i=1:Na
        for j=1:1
            if( B1(1)<Q(i,1) && A3(2)/A3(1)*(Q(i,1)-B3(1))+B3(2) >= Q(i,2) && A4(2)/A4(1)*(Q(i,1)-B4(1))+B4(2) > Q(i,2))
                Mq(mode,i) = Mq(mode,i)+sqrt(1/2/M/Freqs(mode,i))/sqrt(Na)*(-V_ei(Q(i,1),Q(i,2),Q(i,3)))*...
                ( [Q(i,1)-b1(1) Q(i,2)-b1(2) Q(i,3)] * transpose([Eigvecs(mode,i,j,1) Eigvecs(mode,i,j,2) Eigvecs(mode,i,j,3)]));
            elseif(A2(2)/A2(1)*(Q(i,1)-B2(1))+B2(2) < Q(i,2) && A3(2)/A3(1)*(Q(i,1)-B3(1))+B3(2) < Q(i,2) && A5(2)/A5(1)*(Q(i,1)-B5(1))+B5(2) < Q(i,2))
                Mq(mode,i) = Mq(mode,i)+sqrt(1/2/M/Freqs(mode,i))/sqrt(Na)*(-V_ei(Q(i,1),Q(i,2),Q(i,3)))*...
                ( [Q(i,1)-b2(1) Q(i,2)-b2(2) Q(i,3)] * transpose([Eigvecs(mode,i,j,1) Eigvecs(mode,i,j,2) Eigvecs(mode,i,j,3)]));
            elseif( A4(2)/A4(1)*(Q(i,1)-B4(1))+B4(2) <= Q(i,2) && A5(2)/A5(1)*(Q(i,1)-B5(1))+B5(2) > Q(i,2) )
                Mq(mode,i) = Mq(mode,i)+sqrt(1/2/M/Freqs(mode,i))/sqrt(Na)*(-V_ei(Q(i,1),Q(i,2),Q(i,3)))*...
                ( [Q(i,1)-b1(1)-b2(1) Q(i,2)-b1(2)-b2(2) Q(i,3)] * transpose([Eigvecs(mode,i,j,1) Eigvecs(mode,i,j,2) Eigvecs(mode,i,j,3)]));
            else
                Mq(mode,i) = Mq(mode,i)+sqrt(1/2/M/Freqs(mode,i))/sqrt(Na)*(-V_ei(Q(i,1),Q(i,2),Q(i,3)))*...
                ( [Q(i,1) Q(i,2) Q(i,3)] * transpose([Eigvecs(mode,i,j,1) Eigvecs(mode,i,j,2) Eigvecs(mode,i,j,3)]));
            end
        end
    end
end

%%  Cq 함수 개형 보기
Cq=zeros(1,Na);
for i=1:Na
    %V_ei(Q(i,1),Q(i,2),Q(i,3)
    if( B1(1)<Q(i,1) && A3(2)/A3(1)*(Q(i,1)-B3(1))+B3(2) >= Q(i,2) && A4(2)/A4(1)*(Q(i,1)-B4(1))+B4(2) > Q(i,2))
        plot3(Q(i,1)-b1(1),Q(i,2)-b1(2),abs(Mq(5,i))^2/abs(Freqs(4,1)-v*(Q(i,2)-b1(2)))^2,'.','Color','black')
        Cq(i)=abs(Mq(5,i))^2/abs(Freqs(4,1)-v*(Q(i,2)-b1(2)))^2;
        hold on
    elseif(A2(2)/A2(1)*(Q(i,1)-B2(1))+B2(2) < Q(i,2) && A3(2)/A3(1)*(Q(i,1)-B3(1))+B3(2) < Q(i,2) && A5(2)/A5(1)*(Q(i,1)-B5(1))+B5(2) < Q(i,2))
        plot3(Q(i,1)-b2(1),Q(i,2)-b2(2),abs(Mq(5,i))^2/abs(Freqs(4,1)-v*(Q(i,2)-b2(2)))^2,'.','Color','black')
        Cq(i)=abs(Mq(5,i))^2/abs(Freqs(4,1)-v*(Q(i,2)-b2(2)))^2;
        hold on
    elseif( A4(2)/A4(1)*(Q(i,1)-B4(1))+B4(2) <= Q(i,2) && A5(2)/A5(1)*(Q(i,1)-B5(1))+B5(2) > Q(i,2) )
        plot3(Q(i,1)-b1(1)-b2(1),Q(i,2)-b1(2)-b2(2),abs(Mq(5,i))^2/abs(Freqs(4,1)-v*(Q(i,2)-b1(2)-b2(2)))^2,'.','Color','black')
        Cq(i)=abs(Mq(5,i))^2/abs(Freqs(4,1)-v*(Q(i,2)-b1(2)-b2(2)))^2;
        hold on
    else
        plot3(Q(i,1),Q(i,2),abs(Mq(5,i))^2/abs(Freqs(4,1)-v*Q(i,2))^2,'.','Color','black')
        Cq(i)=abs(Mq(5,i))^2/abs(Freqs(4,1)-v*Q(i,2))^2;
        hold on
    end
end
xlabel('Qx')
ylabel('Qy')
title('Cq')
%% Fourier transform
R = 0.4*eV; % Half-Range

N=4000;
Nt=5*N;
t=linspace(0,1000000,Nt);
W=linspace(-R,R,Nt);

Nq=zeros(6,Na);
for mode=6:6
    for q=1:Na
        Nq(mode,q)=1/(exp(Freqs(mode,q)/kb/T)-1);
    end
end


fqy=zeros(Na,1);
for i=1:Na
    if( B1(1)<Q(i,1) && A3(2)/A3(1)*(Q(i,1)-B3(1))+B3(2) >= Q(i,2) && A4(2)/A4(1)*(Q(i,1)-B4(1))+B4(2) > Q(i,2))
        fqy(i)=v*(Q(i,2)-b1(2));
    elseif(A2(2)/A2(1)*(Q(i,1)-B2(1))+B2(2) < Q(i,2) && A3(2)/A3(1)*(Q(i,1)-B3(1))+B3(2) < Q(i,2) && A5(2)/A5(1)*(Q(i,1)-B5(1))+B5(2) < Q(i,2))
        fqy(i)=v*(Q(i,2)-b2(2));
    elseif(A4(2)/A4(1)*(Q(i,1)-B4(1))+B4(2) <= Q(i,2) && A5(2)/A5(1)*(Q(i,1)-B5(1))+B5(2) > Q(i,2) )
        fqy(i)=v*(Q(i,2)-b1(2)-b2(2));
    else
        fqy(i)=v*(Q(i,2));
    end
end

fqy=fqy';

Phif = zeros(1,Nt);
for k=1:(Nt)
    Phif(k) = 0;
    for mode=5:5
        for i=1:Na
            Phif(k) = Phif(k) + abs(Mq(mode,i))^2*(-t(k)*Freqs(mode,i)*1i/(Freqs(mode,i)^2-fqy(i)^2) ...
            + (Nq(mode,i)+1)*(1-exp(-(Freqs(mode,i)-fqy(i))*t(k)*1i))/(Freqs(mode,i)-fqy(i))^2 ...
            + Nq(mode,i)*(1-exp((Freqs(mode,i)+fqy(i))*t(k)*1i))/(Freqs(mode,i)+fqy(i))^2);        
        end
    end
end
% 
G=zeros(1,Nt);

for i=1:Nt
    G(i)=fourier_transform(W(i),Nt,t,Phif,N);
end

Omega = 0;
for q=1:Na
    Omega = Omega+abs(Mq(5,q))^2*Freqs(5,q)/(Freqs(5,q)^2-fqy(q)^2)/eV;
end
shift=(v*0.03663/eV/0.0126-floor(v*0.03663/eV/0.0126))*0.0126
shift=(Omega/0.0126-floor(Omega/0.0126))*0.0126

toc




