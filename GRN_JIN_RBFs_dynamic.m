clc
clear all
%% ��ʼ������
n = 30; %�����˸���
n_tar = 3;
m = 500; %��������
d = 1;
r = 1.5;
ave=0.6283;
%% ��ʼ���ṹ��
z=zeros(n,4,m);
gradient=zeros(n,2,m);
datax=zeros(n,n);
datay=zeros(n,n);
odatax=zeros(n,1);
odatay=zeros(n,1);
wdatax=zeros(n,1);
wdatay=zeros(n,1);
target_position=zeros(n_tar,2,m);
I = zeros(n_tar,2,m);
Total_error=zeros(1,m);
%% ������λ�ó�ʼ��

%% ����һ������ĳ�ʼ����
% z(:,:,1)=8*rand(n,4)-4;
% z(:,:,1)=split_29;
% z(:,:,1)=initial_30_2;
% z(:,:,1)=end2;

%% ����Բ������ĵ���Ϊ��ʼ����
alpha = 2*pi*rand(n,1);
R = 2 + 2*rand(n,1);
xx = R.*cos(alpha);
yy = R.*sin(alpha);
B = zeros(n,2);
% plot(xx,yy,'*');
z(:,:,1) = [xx,yy,B];

%% ��������˶���ʼ��
 target_position(:,:,1) = 4*rand(n_tar,2)-2;
 I (:,:,1)= target_position(:,:,1);
 
%% �ֽ׶��˶���ʼ��
%  target_position(:,:,1) = [0 0;-0.5 0.5;-0.5 -0.5];
%  I (:,:,1)= target_position(:,:,1);


%% �����ʼ���ֵ
for ii=1:n
 val_error = z(ii,1,1).^2+z(ii,2,1).^2-4;   %%Ŀ��ΪԲ�β���
  Error(ii,1) = val_error^2;
end
  Total_error(1,1) = sum(Error,1);


%% ���ѻ�������
% target=[0 0];    %Ŀ����ʼλ��
% plot(0,0,'bp','markerfacecolor','b','markersize',20)
% hold on

% x=0.32;
% target=[2+x 0;-2-x 0;0 2+x;0 -2-x];     %Ŀ����ͼ�뿪���������������
% plot(target(:,1),target(:,2),'bd','markerfacecolor','b','markersize',10)
% hold on
%%  ��ʽ��������
% plot(12,0,'rp','markerfacecolor','r','markersize',10); %����λ��Ϊ��10,11��
% hold on
% plot(8,0,'kp','markersize',10); %�ϰ����ʼλ��Ϊ��10,11��
% hold on
% % %  hold on
% % hold on
%  rectangle('Position',[7.9,2,0.2,4],'facecolor','b');    %��������ʽtunnel
%  hold on
%  rectangle('Position',[7.9,-6,0.2,4],'facecolor','b'); 
%   hold on
%   line([8,12],[0,0]);   
%   hold on

%% ���滷������1
 
% x=0:0.01:20;             
% y=-6*sin(0.1*pi*x)+16;
% plot(x,y,'b');
% hold on  
% 
% line([20,24],[16,16]);         %���ѹ���
% line([24,30],[16,19]);
% line([24,30],[16,13]);
% hold on
% 
% plot([30 30],[19 13],'b.','markerfacecolor','b','markersize',20); %�����ʼλ��Ϊ��10��16���ͣ�11,16��
% hold on
% 
% 
% alpha=0:0.01:2*pi;
% RC=2;
% xc=RC*cos(alpha)+20;
% yc=RC*sin(alpha)+24;
% plot(xc,yc,'-b');
% hold on
%  plot(20,24,'r^','markerfacecolor','r','markersize',10); %�ϰ����ʼλ��Ϊ��20,24��
%  hold on 
%  
% %  rectangle('Position',[30,8,1,16],'facecolor','b'); 
% %  rectangle('Position',[30,8,1,16],'facecolor','b'); 
%  
%  
% rectangle('Position',[30,15.5,5,1],'facecolor','b');      %ǽ������Ϊˮƽ�ģ��̼����ѣ�28,15��
% % 
% 
%% ���滷������2
% x=0:0.01:10;             % ��һ�׶λ�������
% y=-6*sin(0.1*pi*x)+16;
% plot(x,y,'b');
% hold on 
% 
% line([10,36],[16,16]);         %ֱ�������
% line([36,36],[16,20]);
% line([36,36],[16,12]);
% line([36 40],[16 16]);
% hold on
% 
% plot([27.5 28.5 29.5],[16 16 16],'b.','markerfacecolor','b','markersize',20); %�����ʼλ��Ϊ��10��16���ͣ�11,16��
% hold on
% 
% 
% % alpha=0:0.01:2*pi;
% % RC=2;
% % xc=RC*cos(alpha)+10;
% % yc=RC*sin(alpha)+11;
% % plot(xc,yc,'-b');
% % hold on
% %  plot(10,11,'r^','markerfacecolor','r','markersize',10); %�ϰ����ʼλ��Ϊ��10,11��
% %  hold on 
%  
%  rectangle('Position',[16,13,8,1],'facecolor','b');    %�ڶ�����tunnel
%  rectangle('Position',[16,18,8,1],'facecolor','b'); 
%  
%  rectangle('Position',[4,12,2,1],'facecolor','b');    %��һ����tunnel
%  rectangle('Position',[4,7,2,1],'facecolor','b'); 
%  
%  rectangle('Position',[28,6,1,8],'facecolor','b');    %��������ʽtunnel
%  rectangle('Position',[28,18,1,8],'facecolor','b'); 
%  
% % rectangle('Position',[32,11,5,1],'facecolor','b');      %��ǽ������  
% rectangle('Position',[40,12,1,8],'facecolor','b');      %


%% ���滷������3

% x=0:0.01:10;             % ��һ�׶λ�������
% y=-6*sin(0.1*pi*x)+16;
% plot(x,y,'b');
% hold on 
% 
% line([10,36],[16,16]);         %ֱ�������
% % line([36,36],[16,20]);
% % line([36,36],[16,12]);
% % line([36 40],[16 16]);
% hold on
% 
% % %����λ��
% % plot([27.5 28.5 29.5],[16 16 16],'b.','markerfacecolor','b','markersize',20); 
% % hold on
% 
% %����λ��
% plot(34,16,'b.','markerfacecolor','b','markersize',20); 
% hold on
% 
% % z(:,:,1)=robot_32;  % ��һʱ�̻�����λ��
% 
% alpha=0:0.01:2*pi;
% RC=r;
% xc=RC*cos(alpha)+10;
% yc=RC*sin(alpha)+10;
% plot(xc,yc,'-b');
% hold on
% plot(10,10,'r^','markerfacecolor','r','markersize',10);  %�ϰ����ʼλ��Ϊ��10,11��
% hold on 
% 
% 
% 
%  rectangle('Position',[16,13,8,1],'facecolor','b');    %�ڶ�����tunnel
%  rectangle('Position',[16,18,8,1],'facecolor','b'); 
%  
% %  rectangle('Position',[4,12,2,1],'facecolor','b');    %��һ����tunnel
% %  rectangle('Position',[4,7,2,1],'facecolor','b'); 
%  
%  rectangle('Position',[28,6,1,8],'facecolor','b');    %��������ʽtunnel
%  rectangle('Position',[28,18,1,8],'facecolor','b'); 
%  
% % % rectangle('Position',[32,11,5,1],'facecolor','b');      %��ǽ������  
% % rectangle('Position',[40,12,1,8],'facecolor','b');      %


for jj=2:1:m
%% ÿһ��ѭ����ʼǰ�ȼ����������
% d = 0.34 + (1.02-0.34)*(jj/m);   % ���������������
 
%% ����Ŀ��Ķ�̬����  
%  t_trapping=0.01*(jj-1);  %% ������˶��ٶ�Ϊ1m/s   ÿһ��jjѭ����0.2s��ʵ�����ݼ��
%  I = [4*t_trapping 2*sin(t_trapping);3*t_trapping -t_trapping+(0.5*rand()-0.25);3*t_trapping t_trapping+(0.5*rand()-0.25)];
%  target_position(:,:,jj-1)=I;
% [zix,ziy,f]=myfun_zixziy(I);

%% ���Ŀ�겼���˶�����
I(:,:,jj) = I(:,:,jj-1) + (0.2*rand(n_tar,2)-0.1);
target_position(:,:,jj) = I(:,:,jj);
[zix,ziy,f]=myfun_zixziy_any_I(I(:,:,jj-1));


%% ���Ŀ��ֽ׶��˶�����
% if jj <= 100
%     tar1x = 0.01 * (jj-1); tar1y = 0;
%     tar2x = -0.5 + 0.01 * (jj-1); tar2y = 0.5;
%     tar3x = -0.5 + 0.01 * (jj-1); tar3y = -0.5;
% elseif jj<= 200
%     tar1x = 0.01 * (jj-1); tar1y = 0;
%     tar2x = -0.5 + 0.01 * (jj-1); tar2y = 0.5 + 0.03 * (jj-101);
%     tar3x = -0.5 + 0.01 *(jj-1); tar3y = -0.5  - 0.03 * (jj-101);
% elseif jj<=300
%     tar1x = 0.01 * (jj-1); tar1y = 0;
%     tar2x = -0.5 + 0.01 * (jj-1); tar2y = 3.5;
%     tar3x = -0.5 + 0.01 * (jj-1); tar3y = -3.5;
% elseif jj<=400
%     tar1x = 0.01 * (jj-1); tar1y = 0;
%     tar2x = -0.5 + 0.01 * (jj-1); tar2y = 3.5 - 0.03 * (jj-301);
%     tar3x = -0.5 + 0.01 *(jj-1); tar3y = -3.5 + 0.03 * (jj-301);
% else
%     tar1x = 0.01 * (jj-1); tar1y = 0;
%     tar2x = -0.5 + 0.01 * (jj-1); tar2y = 0.5;
%     tar3x = -0.5 + 0.01 * (jj-1); tar3y = -0.5; 
% end
% I(:,:,jj) = [tar1x tar1y;tar2x tar2y;tar3x tar3y];
% target_position(:,:,jj) = I(:,:,jj);
% [zix,ziy,f]=myfun_zixziy_any_I(I(:,:,jj-1));
%% ��������Ϊ��ֱ������������������
%  t_trapping=0.1*(jj-1);  %% ������˶��ٶ�Ϊ1m/s   ÿһ��jjѭ����0.2s��ʵ�����ݼ��
%  I = [0 t_trapping;0 -t_trapping;t_trapping 0];
%  target_position(:,:,jj-1)=I;
% [zix,ziy,f]=myfun_zixziy(I);
%% ���zix ��ziy�Լ�f ͨ������
%  t_trapping=4+0.0*(jj-1);  %% ������˶��ٶ�Ϊ0.1m/s   ÿһ��jjѭ����0.2s��ʵ�����ݼ��
%  I = [t_trapping 0];
%  target_position(:,:,jj-1)=I;
%  [zix,ziy,f]=myfun_zixziy_piece(I,jj);

%%  ����ѭ��
for ii=1:n
    for k=1:n

          if ii==k    %�ų��Լ���Լ�����ɷ�ĸΪ0
              datax(ii,k)=0;
              datay(ii,k)=0;
          else
              
          % �ھӶ���Ϊƽ���������С��d�ĸ��壺
          dij=sqrt((z(k,1,jj-1)-z(ii,1,jj-1))^2+(z(k,2,jj-1)-z(ii,2,jj-1))^2);  
          if dij<d
%              
            valx=(z(ii,1,jj-1)-z(k,1,jj-1))/dij^1;
            valy=(z(ii,2,jj-1)-z(k,2,jj-1))/dij^1;
%             
%             valx=2*(z(ii,1,jj-1)-z(k,1,jj-1))/dij*(dij^2-d0)/dij^3;
%             valy=2*(z(ii,2,jj-1)-z(k,2,jj-1))/dij*(dij^2-d0)/dij^3;
            
%              valx=-2*(z(ii,1,jj-1)-z(k,1,jj-1))/dij^4;
%              valy=-2*(z(ii,2,jj-1)-z(k,2,jj-1))/dij^4;
% %             
      
          
            datax(ii,k)=valx;
            datay(ii,k)=valy;
          else
            datax(ii,k)=0;
            datay(ii,k)=0;
          end
          end
    
    end
        dix=sum(datax,2);
        diy=sum(datay,2);
        
        
         %����ִ���������ϰ������µ���-���꣨10,10��--------------------------------------
           dio=sqrt((z(ii,1,jj-1)-10)^2+(z(ii,2,jj-1)-10)^2);
        if dio<r+d;
            ovalx=(z(ii,1,jj-1)-10)/dio;
            ovaly=(z(ii,2,jj-1)-10)/dio;
            odatax(ii,1)=ovalx;
            odatay(ii,1)=ovaly;
        else
            odatax(ii,1)=0;
            odatay(ii,1)=0;
        end
        odix=sum(odatax,2);
        odiy=sum(odatay,2);
       
        %����ִ���������ϰ������µ���---------------------------------------
%            dio1=sqrt((z(ii,1,jj-1)-5)^2+(z(ii,2,jj-1)-7.5)^2);
%            dio2=sqrt((z(ii,1,jj-1)-5)^2+(z(ii,2,jj-1)-12.5)^2);
%         if dio1<r+d;
%             ovalx=(z(ii,1,jj-1)-5)/dio1;
%             ovaly=(z(ii,2,jj-1)-7.5)/dio1;
%             odatax(ii,1)=ovalx;
%             odatay(ii,1)=ovaly;
%             
%         elseif dio2<r+d;
%             ovalx=(z(ii,1,jj-1)-5)/dio2;
%             ovaly=(z(ii,2,jj-1)-12.5)/dio2;
%             odatax(ii,2)=ovalx;
%             odatay(ii,2)=ovaly;
%         else
%             odatax(ii,1)=0;
%             odatay(ii,1)=0;
%             odatax(ii,2)=0;
%             odatay(ii,2)=0;
%         end
%         odix=sum(odatax,2);
%         odiy=sum(odatay,2);
        
         %����ִ������ǽ�����µ���-----------------------------------------
           diw=sqrt((z(ii,1,jj-1)-40)^2);
        if diw<d;
            wvalx=(z(ii,1,jj-1)-40)/diw;
            wvaly=0;
            wdatax(ii,1)=wvalx;
            wdatay(ii,1)=wvaly;
        else
            wdatax(ii,1)=0;
            wdatay(ii,1)=0;
        end
        wdix=sum(wdatax,2);
        wdiy=sum(wdatay,2);
      

%% �����������ģ������
     x0 = [z(ii,1,jj-1),z(ii,2,jj-1),z(ii,3,jj-1),z(ii,4,jj-1)];
%      t0=0.001*(jj-1):0.005:0.001*jj;
%       t0=0.06*(jj-1):0.03:0.06*jj;
% %      t0=0.06*(jj);
        t0=0.05*(jj-1):0.05:0.05*jj;   %�����ٶ�̫�� ��ô��������������Ŀ��ͼ���ϣ�
%           t0=0.0:0.01:0.2;

%%  ִ����������ͼ����

% ��̬����� Ҫ����zix��ziy�����, ע����myfun_GRN_trapping �м���zix �� ziy����
[t,x]=ode45('myfun_GRN_trapping_dynamic',t0,x0,[],dix(ii),diy(ii),wdix(ii),wdiy(ii),zix,ziy);
  
% ��̬������ ����Ҫ���zix ziy
% [t,x]=ode45('myfun_GRN_trapping',t0,x0,[],dix(ii),diy(ii),wdix(ii),wdiy(ii));
%   
     
     %ִ�л�Բ����
%       xx=34;yy=16;R=4;
%      [t,x]=ode45('myfun_circle_ii_xxyy',t0,x0,[],dix(ii),diy(ii),xx,yy,R,odix(ii),odiy(ii),wdix(ii),wdiy(ii));

    val1=x(end,1); 
    val2=x(end,2);
    val3=x(end,3);
    val4=x(end,4);
%     
%     syms p q 
%     val5=subs(zix,{p,q},{val1+0*val3,val2+0*val4});
%     val6=subs(ziy,{p,q},{val1+0*val3,val2+0*val4});
%     
%     gradient(ii,1,jj)=val5;
%     gradient(ii,2,jj)=val6;
%     
%    val_grad=sqrt(gradient(ii,1,jj)^2+gradient(ii,2,jj)^2);
%    grad(ii,jj)=val_grad;
%    Total_grad=sum(grad,1);
    
    
    z(ii,1,jj)=val1;
    z(ii,2,jj)=val2;
    z(ii,3,jj)=val3;
    z(ii,4,jj)=val4;
    
 %% �����Ƿ���Ŀ��ͼ����
 
  
  %% ����fֻ���ھ�̬Χ��ʱ��Բ���Ŀ��ͼ�ε�����
   syms p q
% ��̬��������2
% f = log((q + 3)^2 + p^2)*((3717804848893897*(q + 3)^2)/36028797018963968 + (3717804848893897*p^2)/36028797018963968) - log((q - 2)^2 + p^2)*((3678916375312215*(q - 2)^2)/9007199254740992 + (3678916375312215*p^2)/9007199254740992) + log((q + 2)^2 + p^2)*((5582798126936243*(q + 2)^2)/18014398509481984 + (5582798126936243*p^2)/18014398509481984) - log((p + 2^(1/2))^2 + (q + 2^(1/2))^2)*((393132877987855*(p + 2^(1/2))^2)/1125899906842624 + (393132877987855*(q + 2^(1/2))^2)/1125899906842624) + log((p - 1)^2 + (q - 1)^2)*((1200804823190509*(p - 1)^2)/9007199254740992 + (1200804823190509*(q - 1)^2)/9007199254740992) + log((p + 1)^2 + (q - 1)^2)*((1200804823190511*(p + 1)^2)/9007199254740992 + (1200804823190511*(q - 1)^2)/9007199254740992) + log((p - 3844062731816691/2251799813685248)^2 + (q - 3844062731816691/2251799813685248)^2)*((5526246206831461*(p - 3844062731816691/2251799813685248)^2)/36028797018963968 + (5526246206831461*(q - 3844062731816691/2251799813685248)^2)/36028797018963968) - log((q + 2^(1/2))^2 + (p - 2^(1/2))^2)*((3145063023902839*(q + 2^(1/2))^2)/9007199254740992 + (3145063023902839*(p - 2^(1/2))^2)/9007199254740992) + log((p + 3844062731816691/2251799813685248)^2 + (q - 3844062731816691/2251799813685248)^2)*((2763123103415731*(p + 3844062731816691/2251799813685248)^2)/18014398509481984 + (2763123103415731*(q - 3844062731816691/2251799813685248)^2)/18014398509481984);
% ��̬��������3
% f = log((p - 8/5)^2 + (q - 9/5)^2)*((4450192366097453*(p - 8/5)^2)/72057594037927936 + (4450192366097453*(q - 9/5)^2)/72057594037927936) - log((q + 4)^2 + p^2)*((6748617341728975*(q + 4)^2)/72057594037927936 + (6748617341728975*p^2)/72057594037927936) + log((p - 1)^2 + (q - 1)^2)*((3060517184830343*(p - 1)^2)/9007199254740992 + (3060517184830343*(q - 1)^2)/9007199254740992) + log((p + 8/5)^2 + (q - 9/5)^2)*((1112548091524361*(p + 8/5)^2)/18014398509481984 + (1112548091524361*(q - 9/5)^2)/18014398509481984) + log((q + 3)^2 + p^2)*((4525188994582225*(q + 3)^2)/9007199254740992 + (4525188994582225*p^2)/9007199254740992) + log((p + 1)^2 + (q - 1)^2)*((6121034369660677*(p + 1)^2)/18014398509481984 + (6121034369660677*(q - 1)^2)/18014398509481984) - log((q - 5/3)^2 + p^2)*((5419212569722201*(q - 5/3)^2)/9007199254740992 + (5419212569722201*p^2)/9007199254740992) - log((p + 6/5)^2 + (q + 29/15)^2)*((2982712687498879*(p + 6/5)^2)/9007199254740992 + (2982712687498879*(q + 29/15)^2)/9007199254740992) - log((p - 6/5)^2 + (q + 29/15)^2)*((2982712687498881*(p - 6/5)^2)/9007199254740992 + (2982712687498881*(q + 29/15)^2)/9007199254740992);
% % ��̬��������4
% f = log((p - 3)^2 + (q + 2)^2)*((5754069591571453*(p - 3)^2)/18014398509481984 + (5754069591571453*(q + 2)^2)/18014398509481984) - log((p - 2509704184906261/4503599627370496)^2 + (q + (6*13^(1/2))/13)^2)*((8722248605890491*(p - 2509704184906261/4503599627370496)^2)/36028797018963968 + (8722248605890491*(q + (6*13^(1/2))/13)^2)/36028797018963968) - log((p + 2)^2 + q^2)*((775656146388575*(p + 2)^2)/9007199254740992 + (775656146388575*q^2)/9007199254740992) - log((p - 11/3)^2 + q^2)*((5059529695416721*(p - 11/3)^2)/18014398509481984 + (5059529695416721*q^2)/18014398509481984) - log((p - 2509704184906261/4503599627370496)^2 + (q - (6*13^(1/2))/13)^2)*((8722248605890481*(p - 2509704184906261/4503599627370496)^2)/36028797018963968 + (8722248605890481*(q - (6*13^(1/2))/13)^2)/36028797018963968) + log((p + 1)^2 + q^2)*((2887990520644761*(p + 1)^2)/9007199254740992 + (2887990520644761*q^2)/9007199254740992) - log((q - 797151290642151/281474976710656)^2 + (p - 8004473239566885/2251799813685248)^2)*((2663551530295575*(q - 797151290642151/281474976710656)^2)/36028797018963968 + (2663551530295575*(p - 8004473239566885/2251799813685248)^2)/36028797018963968) - log((q + 797151290642151/281474976710656)^2 + (p - 8004473239566885/2251799813685248)^2)*((2663551530295607*(q + 797151290642151/281474976710656)^2)/36028797018963968 + (2663551530295607*(p - 8004473239566885/2251799813685248)^2)/36028797018963968) + log((p - 3)^2 + (q - 2)^2)*((2877034795785713*(p - 3)^2)/9007199254740992 + (2877034795785713*(q - 2)^2)/9007199254740992);
%% ����ʱ�̵�������
  val_error= subs(f,{p,q},{z(ii,1,jj),z(ii,2,jj)});
  Error(ii,jj) = val_error^2;
  Total_error = sum(Error,1);
  %% Ŀ��Բ�����仯
%   val_error= z(ii,1,jj)^2+z(ii,2,jj)^2-9;
%   Error(ii,jj) = val_error^2;
%   Total_error = sum(Error,1);
  
  
  
end
fprintf('times: %d \n',jj);
end

%% ���ƻ�������ʼ������λ��
   plot(z(:,1,1),z(:,2,1),'og','LineWidth',2);
   hold on 
   plot(z(:,1,end),z(:,2,end),'or','LineWidth',2);
   hold on 
   
%    plot(0,0,'bp','markerfacecolor','b','markersize',20);
%    
%    hold on

%    plot(c(1:2,1),c(1:2,2),'bd','markerfacecolor','b','markersize',10);
%    hold on
   

   xlabel('x(m)');
   ylabel('y(m)');
   hold on 
%    plot([-5 5 0],[0 0 -5],'bd','markerfacecolor','b','markersize',10);

%% ����������   
 axis equal
%  axis([0 40 0 30]);
   axis([-10 15 -10 10]);
%% ���ƹ켣 ������
 aaa=z(:,1,:);
 bbb=z(:,2,:);
 aaa=reshape(aaa,n,m);
 bbb=reshape(bbb,n,m);
 for kk=1:n
  plot(aaa(kk,:),bbb(kk,:),'-g');
  hold on
 end
%% ����Ŀ������λ��
% I= [1*t_trapping sin(t_trapping);0 -t_trapping;0 t_trapping]; % ��һ������

% I= [4*t_trapping 0;3*t_trapping -t_trapping;3*t_trapping t_trapping]; %����ֱ�߷���

% plot(I(:,1,1),I(:,2,1),'b^','markersize',6);
plot(I(:,1,end),I(:,2,end),'b^','markerfacecolor','b','markersize',6);
%% ����Ŀ��켣
ccc = I(:,1,:);
ddd = I(:,2,:);
ccc = reshape (ccc,n_tar,m);
ddd = reshape (ddd,n_tar,m);
 for kkk=1:n_tar
  plot(ccc(kkk,:),ddd(kkk,:),'-b');
  hold on
 end

%% ����  .mat�ļ�
% save('E:\SCI ˧ ����Χ��\��̬Χ����Ŀ��\dynamic2_data3')