%data 불러오기%
filename1 = 'freqs360000.txt';
filename2 = 'eigvecs360000.txt';
delimiterIn = ' ';
Freqs = importdata(filename1,delimiterIn)/27.211385; %Freqs(mode,num_qs) in atomic unit
eigvecs = importdata(filename2,delimiterIn);
%N은 q point의 개수
dim=size(Freqs);

Eigvecs = zeros(6,Na,2,3);
%zeros(mode,num_qs,atom site,xyz)
for mode=1:6
    for q=1:Na
        for j=1:2
            %x
            Eigvecs(mode,q,j,1)=eigvecs(3*(2-1)+1+6*(mode-1),Na*(j-1)+q)-1i*eigvecs(3*(1-1)+1+6*(mode-1),Na*(j-1)+q);
            %y
            Eigvecs(mode,q,j,2)=eigvecs(3*(2-1)+2+6*(mode-1),Na*(j-1)+q)-1i*eigvecs(3*(1-1)+2+6*(mode-1),Na*(j-1)+q);
            %z
            Eigvecs(mode,q,j,3)=eigvecs(3*(2-1)+3+6*(mode-1),Na*(j-1)+q)-1i*eigvecs(3*(1-1)+3+6*(mode-1),Na*(j-1)+q);
            
            %polarization vector를 phase term을 없애고 재정의 ( normalization은 그대로 유지된다 )
%             Eigvecs(mode,q,j,1)=sqrt(Eigvecs(mode,q,j,1)*conj(Eigvecs(mode,q,j,1)));
%             Eigvecs(mode,q,j,2)=sqrt(Eigvecs(mode,q,j,2)*conj(Eigvecs(mode,q,j,2)));
%             Eigvecs(mode,q,j,3)=sqrt(Eigvecs(mode,q,j,3)*conj(Eigvecs(mode,q,j,3)));
        end
    end
end

for mode=4:6
    for i=1:Na
        Freqs(mode,i) = 0.0126/27.211385;
    end
end

%polarization vector plot
% plot3(Eigvecs(2,:,1,1),Eigvecs(2,:,1,2),Eigvecs(2,:,1,3),'o')
% 
% normalization check
% ans=0;
% for q=1:Na
%     for j=1:1
%         ans=ans + abs(Eigvecs(1,q,j,1))^2+abs(Eigvecs(1,q,j,2))^2+abs(Eigvecs(1,q,j,3))^2;
%     end
% end
% ans
