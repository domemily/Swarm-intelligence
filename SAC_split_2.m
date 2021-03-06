clc
clear all
n=20;
timestep=400;
delta=0.05;
gamma=0.2;
d=3;
d0=1.5;
d1=1.5;
d2=36;

alfa=6;
beta=3;
inta=2;
k_1=3; %雌性个体排斥力强度系数
k_2=2; %雄性个体对雄性个体的作用力系数
k_3=40;%雄性个体对雌性个体的作用力系数
z=zeros(n,5,timestep);



%给定一个随机的初始位置
% a=8*rand(n,2)-4;
% b=0*ones(n,2);
% % c=0*ones(n/10,1);
% % e=1*ones(9*n/10,1);
% ee=zeros(1,1);ff=ones(18,1);
% f=[ee;ff;ee];
% % f=0*ones(n,1);
% z(:,:,1)=[a,b,f];

%给定一个固定的初始位置 雌性个体位于对角线
% a=linspace(0,3,4);
% b=0*ones(4,1);c=1*ones(4,1);e=2*ones(4,1);f=3*ones(4,1);g=4*ones(4,1);
% aa=[b,a';c,a';e,a';f,a';g,a']; ee=zeros(1,1);ff=ones(18,1);gg=ones(2,1);
% bb=zeros(20,2);
% cc=[ee;ff;ee];
% hh=[aa,bb,cc];
% z(:,:,1)=hh;

%给定一个固定的初始位置 雌性个体位于右上角
% a=linspace(0,3,4);
% b=0*ones(4,1);c=1*ones(4,1);e=2*ones(4,1);f=3*ones(4,1);g=4*ones(4,1);
% aa=[b,a';c,a';e,a';f,a';g,a']; ee=zeros(1,1);ff=ones(18,1);gg=ones(2,1);
% bb=zeros(20,2);
% cc=[ff;ee;ee];
% hh=[aa,bb,cc];
% z(:,:,1)=hh;

%给定一个固定的初始位置 雌性个体位于对称边界
a=linspace(0,3,4);
b=0*ones(4,1);c=1*ones(4,1);e=2*ones(4,1);f=3*ones(4,1);g=4*ones(4,1);
aa=[b,a';c,a';e,a';f,a';g,a']; ee=zeros(1,1);ff=ones(7,1);gg=ones(1,1);
bb=zeros(20,2);
cc=[ff;gg;ee;gg;gg;ee;ff;gg];
hh=[aa,bb,cc];
z(:,:,1)=hh;

f=1*ones(n,1);
% zz=[ff;ee;gg;ee;ee;ff];  % 雌性突变个体产生
% zz=f;  % 雌性突变个体产生
zz=[gg;ff;ee;gg;gg;ee;ff;ee];  % 雌性突变个体产生


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
    
%       % 设置速度上限
%         if pix(:)>=3       
%             pix(:)=3;
%         elseif pix(:)<=-3
%             pix(:)=-3;
%         else
%         end
%         
%         if piy(:)>=1 
%             piy(:)=1;
%         elseif piy(:)<=-1
%             piy(:)=-1;
%         else
%         end
% %         %-----------------------------------------

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
          if z(ii,5,jj-1)==0   %如果个体是雌性个体执行下面程序
            if dij<=alfa
                            vvalx=z(k,3,jj-1)-z(ii,3,jj-1);
                            vvaly=z(k,4,jj-1)-z(ii,4,jj-1);
                            sdatax_0(ii,k)=vvalx;
                            sdatay_0(ii,k)=vvaly;
                
              if z(k,5,jj-1)==0
                 valx=k_1*(z(ii,1,jj-1)-z(k,1,jj-1))/dij^1;  %位置协同，同性相斥
                 valy=k_1*(z(ii,2,jj-1)-z(k,2,jj-1))/dij^1;

%                    valx=k_1*((1/dij)-(d2/(dij^3)))*(z(k,1,jj-1)-z(ii,1,jj-1))/dij; %近距离排斥远距离吸引
%                    valy=k_1*((1/dij)-(d2/(dij^3)))*(z(k,2,jj-1)-z(ii,2,jj-1))/dij; %由ii指向k   


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
                  if z(k,5,jj-1)==1
%                       valx=(z(ii,1,jj-1)-z(k,1,jj-1))/dij^2;    %仅有排斥作用
%                       valy=(z(ii,2,jj-1)-z(k,2,jj-1))/dij^2;
                      
                      valx=k_2*((1/dij)-(d0/(dij^3)))*(z(k,1,jj-1)-z(ii,1,jj-1))/dij; %对同性近距离排斥远距离吸引
                      valy=k_2*((1/dij)-(d0/(dij^3)))*(z(k,2,jj-1)-z(ii,2,jj-1))/dij; %由ii指向k   
                      
                      datax_11(ii,k)=valx;
                      datay_11(ii,k)=valy;
                  else 
                      valx=k_3*((1/dij)-(d0/(dij^3)))*(z(k,1,jj-1)-z(ii,1,jj-1))/dij; %对异性近距离排斥远距离吸引
                      valy=k_3*((1/dij)-(d0/(dij^3)))*(z(k,2,jj-1)-z(ii,2,jj-1))/dij; %由ii指向k
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
          
          
% 为速度设置上限---------------------------------
        if pix(ii)>=5       
            pix(ii)=5;
        elseif pix(ii)<=-5
            pix(ii)=-5;
        else
        end
        
        if piy(ii)>=2 
            piy(ii)=2;
        elseif piy(ii)<=-2
            piy(ii)=-2;
        else
        end
        %-----------------------------------------
        
          if z(ii,5,1)==0
              uix(ii)=k_1*femalex_00(ii)+femalex_0(ii)-gamma*pix(ii)+1;
              uiy(ii)=k_1*femaley_00(ii)+femaley_0(ii)-gamma*piy(ii);
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
    
    % 设置速度上限
        if pix(ii)>=7       
            pix(ii)=7;
        elseif pix(ii)<=-7
            pix(ii)=-7;
        else
        end
        
        if piy(ii)>=5 
            piy(ii)=5;
        elseif piy(ii)<=-5
            piy(ii)=-5;
        else
        end
%         %-----------------------------------------

    if jj<=200
    z(ii,5,jj)=z(ii,5,1);
    elseif 200< jj <= 400
    z(ii,5,jj)=zz(ii); 
% z(ii,5,jj)=z(ii,5,1);
    else
    end
   end
   fprintf('times: %d \n',jj);
end

%===============绘制初始位置和最终位置图====================================

%    plot(z(2:19,1,1),z(2:19,2,1),'og');
%    hold on 
%    plot(z(1,1,1),z(1,2,1),'*r');
%    hold on 
%    plot(z(20,1,1),z(20,2,1),'*r');
%    hold on 
%    plot(z(1,1,end),z(1,2,end),'*r','markerfacecolor','r');
%    hold on
%    plot(z(20,1,end),z(20,2,end),'*r','markerfacecolor','r');
%    hold on
%    plot(z(2:19,1,end),z(2:19,2,end),'og','markerfacecolor','g');
 %========================================================================= 
 
 
%===============绘制初始位置和最终位置图====================================
figure(1) 
%起始位置
   plot(z(1:8,1,1),z(1:8,2,1),'og');
   hold on 
   plot(z(10:11,1,1),z(10:11,2,1),'og');
   hold on 
   plot(z(13:20,1,1),z(13:20,2,1),'og');
   hold on 
   plot(z(9,1,1),z(9,2,1),'*r');
   hold on 
   plot(z(12,1,1),z(12,2,1),'*r');
   hold on 
%最终位置
   plot(z(9,1,end),z(9,2,end),'*r','markerfacecolor','r');
   hold on
   plot(z(12,1,end),z(12,2,end),'*r','markerfacecolor','r');
   hold on
   plot(z(n,1,end),z(n,2,end),'*r','markerfacecolor','r');
   hold on
   plot(z(1:8,1,end),z(1:8,2,end),'og','markerfacecolor','g');
   hold on 
   plot(z(10:11,1,end),z(10:11,2,end),'og','markerfacecolor','g');
   hold on 
   plot(z(13:19,1,end),z(13:19,2,end),'og','markerfacecolor','g');
   hold on 
   
   %t=200s的中间位置 新的雌性个体产生
   plot(z(9,1,200),z(9,2,200),'*r','markerfacecolor','r');
   hold on
   plot(z(12,1,200),z(12,2,200),'*r','markerfacecolor','r');
   hold on
   plot(z(n,1,200),z(n,2,200),'*r','markerfacecolor','r');
   hold on
   plot(z(1:8,1,200),z(1:8,2,200),'og','markerfacecolor','g');
   hold on 
   plot(z(10:11,1,200),z(10:11,2,200),'og','markerfacecolor','g');
   hold on 
   plot(z(13:19,1,200),z(13:19,2,200),'og','markerfacecolor','g');
   hold on 
 %=========================================================================  
 
 
 axis equal
 xlabel('x/m');
 ylabel('y/m');


% 绘制轨迹曲线=======================================================

aaa=z(:,1,:);
 bbb=z(:,2,:);
 aaa=reshape(aaa,n,timestep);
 bbb=reshape(bbb,n,timestep);
 for kk=1:n
  plot(aaa(kk,:),bbb(kk,:),':b');
  hold on
 end

%=====================================================================    
    

 ddd=z(:,3,:);
 eee=z(:,4,:);
 aaa=reshape(ddd,n,timestep);
 fff=reshape(eee,n,timestep);
 ccc=1:timestep;
 
 
 figure(2)
 for kk=1:n
  plot(ccc(:),ddd(9,:),'-r');
  hold on
    plot(ccc(:),ddd(12,:),'-r');
  hold on
    plot(ccc(:),ddd(20,:),'-r');
  hold on
    plot(ccc(:),ddd(1:8,:),':b');
  hold on
   plot(ccc(:),ddd(9:11,:),':b');
  hold on
   plot(ccc(:),ddd(13:19,:),':b');
  hold on
 end
 xlabel('t/s');
 ylabel('Vx/(m/s)');
 
figure(3)
 for kk=1:n
  plot(ccc(:),eee(9,:),'-r');
  hold on
    plot(ccc(:),eee(12,:),'-r');
  hold on
    plot(ccc(:),eee(20,:),'-r');
  hold on
    plot(ccc(:),eee(1:8,:),':b');
  hold on
   plot(ccc(:),eee(9:11,:),':b');
  hold on
   plot(ccc(:),eee(13:19,:),':b');
  hold on
 end
 xlabel('t/s)');
 ylabel('Vy/(m/s)');
    
%     
%  figure(4)
%  for kk=1:n
%   plot(aaa(9,:),ddd(9,:),'-r');
%   hold on
%     plot(aaa(12,:),ddd(12,:),'-r');
%   hold on
%     plot(aaa(20,:),ddd(20,:),'-r');
%   hold on
%     plot(aaa(1:8,:),ddd(1:8,:),':b');
%   hold on
%    plot(aaa(9:11,:),ddd(9:11,:),':b');
%   hold on
%    plot(aaa(13:19,:),ddd(13:19,:),':b');
%   hold on
%  end
%  xlabel('x/m)');
%  ylabel('Vx/(m/s)');   
%     
%   figure(5)
%  for kk=1:n
%   plot(bbb(9,:),eee(9,:),'-r');
%   hold on
%     plot(bbb(12,:),eee(12,:),'-r');
%   hold on
%     plot(bbb(20,:),eee(20,:),'-r');
%   hold on
%     plot(bbb(1:8,:),eee(1:8,:),':b');
%   hold on
%    plot(aaa(9:11,:),eee(9:11,:),':b');
%   hold on
%    plot(bbb(13:19,:),eee(13:19,:),':b');
%   hold on
%  end
%  xlabel('y/m)');
%  ylabel('Vy/(m/s)');   
%        
    
    
    
    
    
    