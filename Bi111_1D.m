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
%% delta function
R = 0.4*eV; % Half-Range

N=3000;
Nt=5*N;
t=linspace(0,1000000,Nt);
W=linspace(-R,R,Nt);

%% beta function define
M=3.4702009e-25/9.10938291e-31; %SI/a.u.;
wq=0.0126*eV; 

Nq=1/(exp(wq/kb/T)-1);
a=0.4538*nm;
Na=100000; 
L=(Na-1)*a; % Length of the system

rho=M*Na/L;


q=zeros(1,Na);
for i=1:Na
  q(i)=-pi/a+2*pi/L*(i-1); % generation of q grid
end

fq=v*q;
% Mq=(5e-7)*abs(4*pi*12*e^2./q*sqrt((h/2/pi)/2/rho/L/wq));
% Mq=ones(1,Na)*sqrt(wq^2/Na*3.3); %core model을 보고싶다면
Mq=1/(2*pi)^3*4*pi*83./(q.^2+1^2).*abs(q)*sqrt(1/2/M/wq)/sqrt(Na);
% Mq=(2)*abs(q)*sqrt((h/2/pi)/2/rho/L/wq);

Phif = zeros(1,Nt);
for k=1:(Nt)
    Phif(k) = 0;
    for j=1:Na
        Phif(k) = Phif(k) + Mq(j)^2*(-t(k)*wq*1i/(wq^2-fq(j)^2) ...
        + (Nq+1)*(1-exp(-(wq-fq(j))*t(k)*1i))/((wq-fq(j))^2) ...
        + Nq*(1-exp((wq+fq(j))*t(k)*1i))/((wq+fq(j))^2) );
    end
end
% 
G=zeros(1,Nt);

for i=1:Nt
    G(i)=fourier_transform(W(i),Nt,t,Phif,N);
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

%% plot

plot(W/eV,-2*imag(G))
hold on
toc


