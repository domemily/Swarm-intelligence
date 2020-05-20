clc
clear all
n=20;
timestep=1000;
delta=0.05;
gamma=0.2;
d=3;
d0=1.5;


z=zeros(n,4,timestep);

a=6*rand(n,2)-3;
b=0*ones(n,2);
z(:,:,1)=[a,b];


% z(:,:,1)=6*rand(n,4)-3;

datax=zeros(n,n);
datay=zeros(n,n);
ddatax=zeros(n,n);
ddatay=zeros(n,n);
gix=zeros(n,1);
giy=zeros(n,1);
pix=zeros(n,1);
piy=zeros(n,1);


    gix(:)=z(:,1,1);
    giy(:)=z(:,2,1);
    pix(:)=z(:,3,1);
    piy(:)=z(:,4,1);
    
for jj=2:1:timestep
   for ii=1:n
    for k=1:n

          if ii==k    %�ų��Լ���Լ�����ɷ�ĸΪ0
              datax(ii,k)=0;
              datay(ii,k)=0;
          else
              
          % �ھӶ���Ϊƽ���������С��d�ĸ��壺
          dij=sqrt((z(k,1,jj-1)-z(ii,1,jj-1))^2+(z(k,2,jj-1)-z(ii,2,jj-1))^2);
          if dij<=d
            %����λ��Эͬ��
            valx=2*((1/dij)-(d0/(dij^3)))*(z(k,1,jj-1)-z(ii,1,jj-1))/dij; %�������ų�Զ��������
            valy=2*((1/dij)-(d0/(dij^3)))*(z(k,2,jj-1)-z(ii,2,jj-1))/dij; %��iiָ��k
            %�����ٶ�Эͬ��
            vvalx=z(k,3,jj-1)-z(ii,3,jj-1);
            vvaly=z(k,4,jj-1)-z(ii,4,jj-1);

            datax(ii,k)=valx;
            datay(ii,k)=valy;
            
            ddatax(ii,k)=vvalx;
            ddatay(ii,k)=vvaly;
          else
            datax(ii,k)=0;
            datay(ii,k)=0;
            ddatax(ii,k)=0;
            ddatay(ii,k)=0;
          end
          end
    
    end
        dix=sum(datax,2); %λ��Эͬ�� ��������ӡ�
        diy=sum(datay,2); %λ��Эͬ�� 
        
        speedx=sum(ddatax,2); %�ٶ�Эͬ��
        speedy=sum(ddatay,2); %�ٶ�Эͬ��
        
     %Ϊ�ٶ���������---------------------------------
        if pix(ii)>=4       
            pix(ii)=4;
        elseif pix(ii)<=-4
            pix(ii)=-4;
        else
        end
        
        if piy(ii)>=1 
            piy(ii)=1;
        elseif piy(ii)<=-1
            piy(ii)=-1;
        else
        end
        %-----------------------------------------
        uix(ii)=dix(ii)+speedx(ii)-gamma*pix(ii);
        uiy(ii)=diy(ii)+speedy(ii)-gamma*piy(ii);
        

        
    gix(ii)=gix(ii)+pix(ii)*delta;
    giy(ii)=giy(ii)+piy(ii)*delta;
    pix(ii)=pix(ii)+uix(ii)*delta;
    piy(ii)=piy(ii)+uiy(ii)*delta;
    
    z(ii,1,jj)=gix(ii);
    z(ii,2,jj)=giy(ii);
    z(ii,3,jj)=pix(ii);
    z(ii,4,jj)=piy(ii);
   end
   fprintf('times: %d \n',jj);
end

%===============���Ƴ�ʼλ�ú�����λ��ͼ====================================
%    subplot(1,2,1);
   plot(z(:,1,1),z(:,2,1),'*g');
   hold on 
   plot(z(:,1,end),z(:,2,end),'or');
 %=========================================================================  
%    subplot(1,2,2);
%    plot(t,x);
%    legend('x','y','pix','piy');
%    
 axis equal
%  grid on 

% ���ƹ켣����=======================================================

aaa=z(:,1,:);
 bbb=z(:,2,:);
 aaa=reshape(aaa,n,timestep);
 bbb=reshape(bbb,n,timestep);
 for kk=1:n
  plot(aaa(kk,:),bbb(kk,:),':b');
  hold on
 end
%  axis([-2 2 -2 2]);
%   axis([-5 5 -5 5]);
%=====================================================================    
    
    
    
    
    
    
    
    
    
    
    
    
    