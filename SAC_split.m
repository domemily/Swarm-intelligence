clc
clear all
n=20;
timestep=200;
delta=0.05;
gamma=0.2;
d=3;
d0=1.5;
alfa=6;
beta=2;
inta=2;
z=zeros(n,5,timestep);

a=6*rand(n,2)-3;
b=0*ones(n,2);
c=1*ones(n/10,1);
e=0*ones(9*n/10,1);
% f=[c;e];
f=0*ones(n,1);
z(:,:,1)=[a,b,f];


% z(:,:,1)=6*rand(n,4)-3;

datax_00=zeros(n,n);
datay_00=zeros(n,n);
datax_01=zeros(n,n);
datay_01=zeros(n,n);
datax_11=zeros(n,n);
datay_11=zeros(n,n);
datax_10=zeros(n,n);
datay_10=zeros(n,n);

sdatax_0=zeros(n,n);
sdatay_0=zeros(n,n);
sdatax_1=zeros(n,n);
sdatay_1=zeros(n,n);

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

          if ii==k    %排除自减项，自减会造成分母为0
                datax_00(ii,k)=0;
                datay_00(ii,k)=0;
                datax_01(ii,k)=0;
                datay_01(ii,k)=0;
                datax_11(ii,k)=0;
                datay_11(ii,k)=0;
                datax_10(ii,k)=0;
                datay_10(ii,k)=0;
                
          else
              
          % 邻居定义为平面坐标距离小于d的个体：
          dij=sqrt((z(k,1,jj-1)-z(ii,1,jj-1))^2+(z(k,2,jj-1)-z(ii,2,jj-1))^2);
          if z(ii,5,1)==0   %如果个体是雌性个体执行下面程序
            if dij<=alfa
                            vvalx=z(k,3,jj-1)-z(ii,3,jj-1);
                            vvaly=z(k,4,jj-1)-z(ii,4,jj-1);
                            sdatax_0(ii,k)=vvalx;
                            sdatay_0(ii,k)=vvaly;
                
              if z(k,5,1)==0
                 valx=(z(ii,1,jj-1)-z(k,1,jj-1))/dij^2;  %位置协同，同性相斥
                 valy=(z(ii,2,jj-1)-z(k,2,jj-1))/dij^2;
                 

                  datax_00(ii,k)=valx;
                  datay_00(ii,k)=valy;
                  

              else 
                  valx=z(k,1,jj-1)-z(ii,1,jj-1);       
                  valy=z(k,2,jj-1)-z(ii,2,jj-1);
                  datax_01(ii,k)=valx;
                  datay_01(ii,k)=valy;
              end
            else
                datax_00(ii,k)=0;
                datay_00(ii,k)=0;
                datax_01(ii,k)=0;
                datay_01(ii,k)=0;
                
                sdatax_0(ii,k)=0;
                sdatay_0(ii,k)=0;
            end
          else        %如果个体是雄性个体执行下面程序
              if dij<=beta
                            vvalx=z(k,3,jj-1)-z(ii,3,jj-1);
                            vvaly=z(k,4,jj-1)-z(ii,4,jj-1);
                            sdatax_1(ii,k)=vvalx;
                            sdatay_1(ii,k)=vvaly;
                  if z(k,5,1)==1
                      valx=(z(ii,1,jj-1)-z(k,1,jj-1))/dij^2;    %仅有排斥作用
                      valy=(z(ii,2,jj-1)-z(k,2,jj-1))/dij^2;
                      
%                       valx=2*((1/dij)-(d0/(dij^3)))*(z(k,1,jj-1)-z(ii,1,jj-1))/dij; %近距离排斥远距离吸引
%                       valy=2*((1/dij)-(d0/(dij^3)))*(z(k,2,jj-1)-z(ii,2,jj-1))/dij; %由ii指向k   
                      
                      datax_11(ii,k)=valx;
                      datay_11(ii,k)=valy;
                  else 
                      valx=2*((1/dij)-(d0/(dij^3)))*(z(k,1,jj-1)-z(ii,1,jj-1))/dij; %近距离排斥远距离吸引
                      valy=2*((1/dij)-(d0/(dij^3)))*(z(k,2,jj-1)-z(ii,2,jj-1))/dij; %由ii指向k
                      datax_10(ii,k)=valx;
                      datay_10(ii,k)=valy;
                  end
              else
                  datax_11(ii,k)=0;
                  datay_11(ii,k)=0;
                  datax_10(ii,k)=0;
                  datay_10(ii,k)=0;
                  
                  sdatax_1(ii,k)=0;
                  sdatay_1(ii,k)=0;
                  
                  
              end
          end
          end
    end
          femalex_00=sum(datax_00,2);
          femaley_00=sum(datay_00,2);
          femalex_01=sum(datax_01,2);
          femaley_01=sum(datay_01,2);
          
          femalex_0=sum(sdatax_0,2);
          femaley_0=sum(sdatay_0,2);
          
          malex_11=sum(datax_11,2);
          maley_11=sum(datay_11,2);
          malex_10=sum(datax_10,2);
          maley_10=sum(datay_10,2);
             
          malex_1=sum(sdatax_1,2);
          maley_1=sum(sdatay_1,2);
          
          
               %为速度设置上限---------------------------------
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
        
          if z(ii,5,1)==0
              uix(ii)=femalex_00(ii)+inta*rand()*femalex_01(ii)+femalex_0(ii)-gamma*pix(ii)+1;
              uiy(ii)=femaley_00(ii)+inta*rand()*femaley_01(ii)+femaley_0(ii)-gamma*piy(ii);
          else 
              uix(ii)=malex_11(ii)+malex_10(ii)+malex_1(ii)-gamma*pix(ii)+1;
              uiy(ii)=maley_11(ii)+maley_10(ii)+maley_1(ii)-gamma*piy(ii);
          end
          
        
     
    gix(ii)=gix(ii)+pix(ii)*delta;
    giy(ii)=giy(ii)+piy(ii)*delta;
    pix(ii)=pix(ii)+uix(ii)*delta;
    piy(ii)=piy(ii)+uiy(ii)*delta;
    
    z(ii,1,jj)=gix(ii);
    z(ii,2,jj)=giy(ii);
    z(ii,3,jj)=pix(ii);
    z(ii,4,jj)=piy(ii);
    z(ii,5,jj)=z(ii,5,1);
   end
   fprintf('times: %d \n',jj);
end

%===============绘制初始位置和最终位置图====================================
%    subplot(1,2,1);
   plot(z(:,1,1),z(:,2,1),'ok');
   hold on 
   plot(z(1:n/10,1,end),z(1:n/10,2,end),'or');
   hold on
   plot(z(n/10+1:n,1,end),z(n/10+1:n,2,end),'og');
 %=========================================================================  
%    subplot(1,2,2);
%    plot(t,x);
%    legend('x','y','pix','piy');
%    
 axis equal
%  grid on 

% 绘制轨迹曲线=======================================================

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
    
    
    
    
    
    
    
    
    
    
    
    
    