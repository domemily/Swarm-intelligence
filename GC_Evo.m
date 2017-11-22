clc
clear
 rand('seed',15) 
%% initial paremeters
n_prey=2;
n_pred=2;
timestep=1000; 
delta=0.05;

alfa = 1; % for speed match
beta=1;   % for position match
gamma=5.2;  % predation interaction
omg=5; % wall interaction

v_co = 10;  %max speed of chaser
v_to = 10;   % max speed of target
k_chase = -0.001; % the parameter of energy lose of predator during chasing
k_escape = -0.001;  %the parameter of energy lose of prey during escaping
k_nat = -0.002; % natural degradation
r_death = 1; %
delay= 20;
d0=1;   % expected distance

L=20;
R_area=20;  %radius of simulation world
R_wall=0.5; % width of wall
R_nei_prey = 5*r_death;  % neighbor distance for local best of prey
R_nei_pred = 5*r_death; % neighbor distance for local best of pred
U_c=0.0001; % the energy obtaining rate for a neighborless prey
unit_s = 0.1; % the energy obtained by pred from each prey when it is dead
in_e = 0.8; % initial energy
max_e = 1; % max energy
min_e = 0.3; %min energy
noise_pred = 0.1; %
noise_prey = 0.1; %
noise_active = 1; % start again with another random direction
cog = 1; % cognitive weight
soc = 1; % social weight
%% initial structure
z_prey=zeros(n_prey,10,timestep);
z_pred=zeros(n_pred,10,timestep);

pos_datax_22=zeros(n_pred,n_pred);
pos_datay_22=zeros(n_pred,n_pred);
spe_datax_22=zeros(n_pred,n_pred);
spe_datay_22=zeros(n_pred,n_pred);
pos_datax_11=zeros(n_prey,n_prey);
pos_datay_11=zeros(n_prey,n_prey);
spe_datax_11=zeros(n_prey,n_prey);
spe_datay_11=zeros(n_prey,n_prey);
datax_11=zeros(n_prey,n_prey);
datay_11=zeros(n_prey,n_prey);
datax_pred=zeros(n_pred,n_prey);
datay_pred=zeros(n_pred,n_prey);
chase_21=zeros(n_pred,n_prey);
chase_12=zeros(n_prey,n_pred);
f_ct_x=zeros(n_pred,1);
f_ct_y=zeros(n_pred,1);
prey_worse=zeros(n_prey,n_pred);
prey_normal=zeros(n_prey,n_pred);
stop_pred=zeros(n_pred,timestep);
E_ind=zeros(n_pred,n_prey);
gix_2=zeros(n_pred,1);
giy_2=zeros(n_pred,1);
pix_2=zeros(n_pred,1);
piy_2=zeros(n_pred,1);
gix_1=zeros(n_prey,1);
giy_1=zeros(n_prey,1);
pix_1=zeros(n_prey,1);
piy_1=zeros(n_prey,1);
fit_prey=zeros(n_prey,timestep);
sbest_prey=zeros(n_prey,6,timestep);
lbest_prey=zeros(n_prey,6,timestep);
fit_pred=zeros(n_pred,timestep);
sbest_pred=zeros(n_pred,6,timestep);
lbest_pred=zeros(n_pred,6,timestep);
%% initial position
%------------------------square initial postion--------------------
% a=L*rand(n_prey,4)-L/2;b=L*rand(n_pred,4)-L/2;c=in_e*ones(n_prey,1);d=in_e*ones(n_pred,1);e=rand(n_prey,5);f=rand(n_pred,5);
% z_pred(:,:,1)=[b,d,f];
% z_prey(:,:,1)=[a,c,e];

%------------ circle initial position--------------------------------------

ang_prey = 2*pi*rand(n_prey,1);ang_pred = 2*pi*rand(n_pred,1);
R_prey = R_area*rand(n_prey,1);R_pred = R_area*rand(n_pred,1);
xx = R_prey.*cos(ang_prey); xxx = R_pred.*cos(ang_pred);
yy = R_prey.*sin(ang_prey); yyy = R_pred.*sin(ang_pred);
a = [xx,yy]; aa=[xxx,yyy];
 b=zeros(n_prey,2); c = in_e*ones(n_prey,1);d=rand(n_prey,5);
 bb=zeros(n_pred,2); cc = in_e*ones(n_pred,1);dd=rand(n_pred,5);

z_prey(:,:,1)=[a,b,c,d];
z_pred(:,:,1)=[aa,bb,cc,dd];

%% initial positon and veclocity
gix_2(:)=z_pred(:,1,1);
giy_2(:)=z_pred(:,2,1);
pix_2(:)=z_pred(:,3,1);
piy_2(:)=z_pred(:,4,1);
gix_1(:)=z_prey(:,1,1);
giy_1(:)=z_prey(:,2,1);
pix_1(:)=z_prey(:,3,1);
piy_1(:)=z_prey(:,4,1);

%% initial fit, sbest and gbest
fit_prey(:,1)=zeros(n_prey,1); fit_pred(:,1)=zeros(n_pred,1);
sbest_prey(:,6,1)=zeros(n_prey,1); sbest_pred(:,6,1)=zeros(n_pred,1);
sbest_prey(:,1:5,1)=z_prey(:,6:10,1);sbest_pred(:,1:5,1)=z_pred(:,6:10,1);
lbest_prey(:,1:5,1) = z_prey(:,6:10,1);lbest_pred(:,1:5,1) = z_pred(:,6:10,1); % asume 1st agent is the global best
lbest_prey(:,6,1) = zeros(n_prey,1);lbest_pred(:,6,1) = zeros(n_pred,1);

%% main loop

for jj=2:timestep
  
    %% programme of predator
    for ii=1:n_pred
  
        for k=1:n_pred
 
          d_22=sqrt((z_pred(k,1,jj-1)-z_pred(ii,1,jj-1))^2+...
              (z_pred(k,2,jj-1)-z_pred(ii,2,jj-1))^2);
          if d_22<=z_pred(ii,7,jj-1)*R_area          
            if d_22==0
            pos_datax_22(ii,k)=0;
            pos_datay_22(ii,k)=0; 
            else
            valx=2*((1/d_22)-(d0/(d_22^3)))*(z_pred(k,1,jj-1)-z_pred(ii,1,jj-1))/d_22; %
            valy=2*((1/d_22)-(d0/(d_22^3)))*(z_pred(k,2,jj-1)-z_pred(ii,2,jj-1))/d_22;
            pos_datax_22(ii,k)=valx;
            pos_datay_22(ii,k)=valy;
            end
          else
            pos_datax_22(ii,k)=0;
            pos_datay_22(ii,k)=0;
          end
          if d_22<=z_pred(ii,6,jj-1)*R_area 
            vvalx=z_pred(k,3,jj-1);
            vvaly=z_pred(k,4,jj-1);
            spe_datax_22(ii,k)=vvalx;
            spe_datay_22(ii,k)=vvaly;
          else
            spe_datax_22(ii,k)=0;
            spe_datay_22(ii,k)=0;
          end
    
       end
    pos_x_22=sum(pos_datax_22,2);
    pos_y_22=sum(pos_datay_22,2);
    spe_x_22=sum(spe_datax_22,2);
    spe_y_22=sum(spe_datay_22,2);
        
    %% interaction between two groups for predator
        for kk=1:n_prey
            
           if z_prey(kk,5,jj-1)==-1000  % the prey is dead
           chase_21(ii,kk)=1000;  % make no influence 
           datax_pred(ii,kk) = 0; 
           datay_pred(ii,kk) = 0;
           E_ind(ii,kk)=0; % the number of dead
           else
        d_21=sqrt((z_pred(ii,1,jj-1)-z_prey(kk,1,jj-1))^2+...
              (z_pred(ii,2,jj-1)-z_prey(kk,2,jj-1))^2);
          if d_21<= r_death 
              chase_21(ii,kk)=1000;
              E_ind(ii,kk)=1; % the number of dead
              valx = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,1,jj-1)-z_prey(kk,1,jj-1));
              valy = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,2,jj-1)-z_prey(kk,2,jj-1));
              datax_pred(ii,kk)=valx;
              datay_pred(ii,kk)=valy; 
          elseif d_21<= z_pred(ii,9,jj-1)*R_area && d_21 >r_death
              val = d_21;
              chase_21(ii,kk)=val;
              valx = 1/(chase_21(ii,kk))^2*...
                     (z_pred(ii,1,jj-1)-z_prey(kk,1,jj-1));
              valy = 1/(chase_21(ii,kk))^2*...
                     (z_pred(ii,2,jj-1)-z_prey(kk,2,jj-1));
              datax_pred(ii,kk)=valx;
              datay_pred(ii,kk)=valy; 
              E_ind(ii,kk)=0; % the number of dead
          else
              E_ind(ii,kk)=0; % the number of dead
              chase_21(ii,kk)=1000;
              valx = 1/(chase_21(ii,kk))^10;
              valy = 1/(chase_21(ii,kk))^10;
              datax_pred(ii,kk)=valx;
              datay_pred(ii,kk)=valy;
              
          end
          end
        end
        
        
        %% arena
         d_area=sqrt(z_pred(ii,1,jj-1)^2+z_pred(ii,2,jj-1)^2);
         if d_area<R_area
             s=0;
         elseif d_area>=R_area && d_area<=R_area+R_wall
             s=-1/2*sin(pi*(d_area-R_area)/R_wall-pi/2)-1/2;
         elseif d_area>R_area+R_wall
             
             s=-1;
         else
       
         end
         
        f_wall_x(ii)=s*(v_co*z_pred(ii,1,jj-1)/d_area+z_pred(ii,3,jj-1));
        f_wall_y(ii)=s*(v_co*z_pred(ii,2,jj-1)/d_area+z_pred(ii,4,jj-1));
        
          ct_x(ii)=sum(datax_pred(ii,:));
          ct_y(ii)=sum(datay_pred(ii,:));
          f_ct_x(ii)=ct_x(ii)/sqrt((ct_x(ii))^2+(ct_y(ii))^2);  % 
          f_ct_y(ii)=ct_y(ii)/sqrt((ct_x(ii))^2+(ct_y(ii))^2);
              
          uix_2(ii)=beta*pos_x_22(ii)+alfa*spe_x_22(ii)-gamma*f_ct_x(ii)+omg*f_wall_x(ii);
          uiy_2(ii)=beta*pos_y_22(ii)+alfa*spe_y_22(ii)-gamma*f_ct_y(ii)+omg*f_wall_y(ii);
          
          v_c(ii) = z_pred(ii,8,jj-1)*v_co*sigmf(z_pred(ii,5,jj-1),[10 0.5]); % The magnitude of the velocity, a variable
%         v_c(ii) = v_co;
          % speed update
          pix_2(ii)=v_c(ii)*uix_2(ii)/(sqrt(uix_2(ii)^2+uiy_2(ii)^2));
          piy_2(ii)=v_c(ii)*uiy_2(ii)/(sqrt(uix_2(ii)^2+uiy_2(ii)^2));
          % noise
          xita_pred = 2*noise_pred*pi*rand()-noise_pred*pi;
          pix_2(ii)=pix_2(ii)*cos(xita_pred)-piy_2(ii)*sin(xita_pred);
          piy_2(ii)=piy_2(ii)*cos(xita_pred)+pix_2(ii)*sin(xita_pred);

          gix_2(ii)=gix_2(ii)+pix_2(ii)*delta;
          giy_2(ii)=giy_2(ii)+piy_2(ii)*delta;
          
    z_pred(ii,1,jj)=gix_2(ii);
    z_pred(ii,2,jj)=giy_2(ii);
    z_pred(ii,3,jj)=pix_2(ii);
    z_pred(ii,4,jj)=piy_2(ii);
    z_pred(ii,6:10,jj)=z_pred(ii,6:10,jj-1);
%% obtain theta
% dot(ii) = z_pred(ii,3,jj-1)*pix_2(ii)+z_pred(ii,4,jj-1)*piy_2(ii);
% module(ii) = sqrt((z_pred(ii,3,jj-1))^2+(z_pred(ii,4,jj-1))^2)*...
% sqrt((pix_2(ii))^2+(piy_2(ii))^2);
% yuxuan(ii)=dot(ii)/module(ii);            % cos(theta)=a*b/(/a/*/b/)
% theta(ii) = acos(yuxuan(ii));


dot(ii)=pix_2(ii)-z_pred(ii,3,jj-1);
module(ii)=sqrt((pix_2(ii)-z_pred(ii,3,jj-1))^2+(piy_2(ii)-z_pred(ii,4,jj-1))^2);
yuxuan(ii)=dot(ii)/module(ii);            % cos(theta)=a/sart(a^2+b^2)
theta(ii) = acos(yuxuan(ii));

%% energy update
if z_pred(ii,5,jj-1)<min_e 
    z_pred(ii,5,jj)=-1000;  % dead label
elseif z_pred(ii,5,jj-1)>=max_e
    z_pred(ii,5,jj)=max_e;
else
end   


 
%% obtain fit, sbest and gbest of pred

if z_pred(ii,5,jj)==-1000 % obtain fit
    fit_pred(ii,jj)=0;
else 
    fit_pred(ii,jj)=fit_pred(ii,jj-1)+1;
end

if fit_pred(ii,jj)>sbest_pred(ii,6,jj-1) % obtain sbest
    sbest_pred(ii,6,jj)=fit_pred(ii,jj);
    sbest_pred(ii,1:5,jj)=z_pred(ii,6:10,jj);
else 
    sbest_pred(ii,:,jj)=sbest_pred(ii,:,jj-1);
    
end

for k = 1 :n_pred   %obtain local best
   d_22=sqrt((z_pred(k,1,jj-1)-z_pred(ii,1,jj-1))^2+...
              (z_pred(k,2,jj-1)-z_pred(ii,2,jj-1))^2);
          if d_22<R_nei_pred
              val = fit_pred(k,jj);
              fit_pred_neigh(ii,k) = val;
          else
              val = 0;
              fit_pred_neigh(ii,k) = val; 
          end

end

[lbest_pred(ii,6,jj),index]=max(fit_pred_neigh(ii,:));  % locl best

if lbest_pred(ii,6,jj)>lbest_pred(ii,6,jj-1)
    lbest_pred(ii,6,jj)=lbest_pred(ii,6,jj);
    lbest_pred(ii,1:5,jj)=z_pred(index,6:10,jj-1);
else
    lbest_pred(ii,6,jj)=lbest_pred(ii,6,jj-1);
    lbest_pred(ii,:,jj)=lbest_pred(ii,:,jj-1); 
end


% rebirth of predators and energy update

if z_pred(ii,5,jj)==-1000
    
   %---------------particle swarm optimization tendency------
    A = rand();B = rand(); C=A+B;
   z_pred(ii,6,jj) = z_pred(ii,6,jj-1)+...
     cog*(A/C)*(sbest_pred(ii,1,jj)-z_pred(ii,6,jj-1))+...
     soc*(B/C)*(lbest_pred(ii,1,jj)-z_pred(ii,6,jj-1));
   z_pred(ii,7,jj) = z_pred(ii,7,jj-1)+...
     cog*(A/C)*(sbest_pred(ii,2,jj)-z_pred(ii,7,jj-1))+...
     soc*(B/C)*(lbest_pred(ii,2,jj)-z_pred(ii,7,jj-1));
   z_pred(ii,8,jj) = z_pred(ii,8,jj-1)+...
     cog*(A/C)*(sbest_pred(ii,3,jj)-z_pred(ii,8,jj-1))+...
     soc*(B/C)*(lbest_pred(ii,3,jj)-z_pred(ii,8,jj-1));
   z_pred(ii,9,jj) = z_pred(ii,9,jj-1)+...
     cog*(A/C)*(sbest_pred(ii,4,jj)-z_pred(ii,9,jj-1))+...
     soc*(B/C)*(lbest_pred(ii,4,jj)-z_pred(ii,9,jj-1));
   z_pred(ii,10,jj) = z_pred(ii,10,jj-1)+...
     cog*(A/C)*(sbest_pred(ii,5,jj)-z_pred(ii,10,jj-1))+...
     soc*(B/C)*(lbest_pred(ii,5,jj)-z_pred(ii,10,jj-1));
    
 %---------------%average tendency------    
    
%   AVG_pred = (sum(z_pred(:,:,jj-1),1)-z_pred(ii,:,jj-1))/(n_pred-1);
%   z_pred(ii,6:10,jj) = AVG_pred(6:10);    
    
 %----------------------------------------------------------   
    z_pred(ii,5,jj)=0.8;   % the energy level of newbirth
  
    angle = 2*pi*rand();
    R_re = R_area*rand();
    a = R_re*cos(angle);
    b = R_re*sin(angle);
    c=zeros(1,2);
    z_pred(ii,1:4,jj)=[a,b,c];    
    
    gix_2(ii)=z_pred(ii,1,jj);
    giy_2(ii)=z_pred(ii,2,jj);
    pix_2(ii)=z_pred(ii,3,jj);
    piy_2(ii)=z_pred(ii,4,jj);
else 
    
E_add(ii) = sum(E_ind(ii,:))*unit_s;%
z_pred(ii,5,jj) = z_pred(ii,5,jj-1) + E_add(ii) +...
k_chase*((v_c(ii)/v_co)^2+(theta(ii)/pi)^2)+...
k_nat; % 
end



%% handling time 
if jj>3
if z_pred(ii,5,jj-1)>z_pred(ii,5,jj-2) 
    stop_pred(ii,jj)=jj-1;
elseif z_pred(ii,5,jj-1)==z_pred(ii,5,jj-2) 
    stop_pred(ii,jj)=stop_pred(ii,jj-1);
else
    stop_pred(ii,jj)=jj-1-delay;
end
if jj<stop_pred(ii,jj)+delay
    z_pred(ii,:,jj)=z_pred(ii,:,jj-1);
    gix_2(ii) = z_pred(ii,1,jj); % update gix and giy 
    giy_2(ii) = z_pred(ii,2,jj);  % 
    pix_2(ii) = z_pred(ii,3,jj);
    piy_2(ii) = z_pred(ii,4,jj);
    % random moving direction when start again
    xita_active= 2*noise_active*pi*rand()-noise_active*pi;
    pix_2(ii)=pix_2(ii)*cos(xita_active)-piy_2(ii)*sin(xita_active);
    piy_2(ii)=piy_2(ii)*cos(xita_active)+pix_2(ii)*sin(xita_active);
    z_pred(ii,3,jj)=pix_2(ii);
    z_pred(ii,4,jj)=piy_2(ii);
    
end
end   
          
    
    end  
    
%% the stage of prey==========================================================
   for ii=1:n_prey
         
            
      for k=1:n_prey   % interaction among preys
          % 
          d_11=sqrt((z_prey(k,1,jj-1)-z_prey(ii,1,jj-1))^2+...
              (z_prey(k,2,jj-1)-z_prey(ii,2,jj-1))^2);
          if d_11<=z_prey(ii,7,jj-1)*R_area
            %
            if d_11==0
                pos_datax_11(ii,k)=0;
                pos_datay_11(ii,k)=0;
                prey_neibour(ii,k)=1;
            else
            valx=2*((1/d_11)-(d0/(d_11^3)))*(z_prey(k,1,jj-1)-z_prey(ii,1,jj-1))/d_11; %
            valy=2*((1/d_11)-(d0/(d_11^3)))*(z_prey(k,2,jj-1)-z_prey(ii,2,jj-1))/d_11;
            pos_datax_11(ii,k)=valx;
            pos_datay_11(ii,k)=valy;
            prey_neibour(ii,k)=1; % obtain the competition from other preys, neighos come from postin match
            end
          else
            pos_datax_11(ii,k)=0;
            pos_datay_11(ii,k)=0;
            prey_neibour(ii,k)=0; % obtain the competition from other preys
          end
          if d_11<=z_prey(ii,6,jj-1)*R_area %
            vvalx=z_prey(k,3,jj-1);
            vvaly=z_prey(k,4,jj-1);
            spe_datax_11(ii,k)=vvalx;
            spe_datay_11(ii,k)=vvaly;
          else
            spe_datax_11(ii,k)=0;
            spe_datay_11(ii,k)=0;
          end
       
      end
    pos_x_11=sum(pos_datax_11,2);
    pos_y_11=sum(pos_datay_11,2);
    spe_x_11=sum(spe_datax_11,2);
    spe_y_11=sum(spe_datay_11,2);
    Peer_num = sum(prey_neibour,2); % obtain the competition from other preys, neighos come from postin match
   
    %% interaction between two groups for prey       
    for kk=1:n_pred      
        if z_prey(ii,5,jj-1)==-1000
           chase_12(ii,kk)=1000;  % make no influence 
           datax(ii,kk) = 0; 
           datay(ii,kk) = 0;
        else
 
         d_12=sqrt((z_prey(ii,1,jj-1)-z_pred(kk,1,jj-1))^2+...
              (z_prey(ii,2,jj-1)-z_pred(kk,2,jj-1))^2);
          if d_12<=z_prey(ii,9,jj-1)*R_area && d_12>=r_death  % see the predator but alive
              val = d_12;
              chase_12(ii,kk)=val;
              valx = 1/(chase_12(ii,kk)-r_death)^2*...
                     (z_prey(ii,1,jj-1)-z_pred(kk,1,jj-1));
              valy = 1/(chase_12(ii,kk)-r_death)^2*...
                     (z_prey(ii,2,jj-1)-z_pred(kk,2,jj-1));
              datax(ii,kk)=valx;
              datay(ii,kk)=valy;
              prey_normal(ii,kk)=0;
              prey_worse(ii,kk)=1;
          elseif  d_12<r_death   
              z_prey(ii,5,jj)=-1000; % label the dead, the only place to deal with the dead prey
              chase_12(ii,kk)=1000;
              valx = 1/(chase_12(ii,kk)-r_death)^2*...
                     (z_prey(ii,1,jj-1)-z_pred(kk,1,jj-1));
              valy = 1/(chase_12(ii,kk)-r_death)^2*...
                     (z_prey(ii,2,jj-1)-z_pred(kk,2,jj-1));
              datax(ii,kk)=valx;
              datay(ii,kk)=valy;  
              prey_normal(ii,kk)=1;
              prey_worse(ii,kk)=0;
          else % 
              chase_12(ii,kk)=1000;
              valx = 1/(chase_12(ii,kk)-r_death)^2;
              valy = 1/(chase_12(ii,kk)-r_death)^2;
              datax(ii,kk)=valx;
              datay(ii,kk)=valy; 
              prey_normal(ii,kk)=1;
              prey_worse(ii,kk)=0;
          end
        end
    end
    worse(ii)=any(prey_worse(ii,:)==1); % logic of run
    normal(ii)=all(prey_normal(ii,:)==1); % logic of eating
    prey_stop(ii)=1-normal(ii); %logic of prey moving 
    %% arean constraints
    %% arena
         d_area=sqrt(z_prey(ii,1,jj-1)^2+z_prey(ii,2,jj-1)^2);
         if d_area<R_area
             s=0;
         elseif d_area>=R_area && d_area<=R_area+R_wall
             s=-1/2*sin(pi*(d_area-R_area)/R_wall-pi/2)-1/2;
         elseif d_area>R_area+R_wall
             
             s=-1;
         else
       
         end
         
        f_wall_x(ii)=s*(v_to*z_prey(ii,1,jj-1)/d_area+z_prey(ii,3,jj-1));
        f_wall_y(ii)=s*(v_to*z_prey(ii,2,jj-1)/d_area+z_prey(ii,4,jj-1));
    

              ctx(ii)=sum(datax(ii,:));
              cty(ii)=sum(datay(ii,:));
              f_ct_x(ii)=ctx(ii)/sqrt((ctx(ii))^2+(cty(ii))^2);  % 
              f_ct_y(ii)=cty(ii)/sqrt((ctx(ii))^2+(cty(ii))^2);
           
   %% 
          uix_1(ii)=beta*pos_x_11(ii)+alfa*spe_x_11(ii)+gamma*f_ct_x(ii)+omg*f_wall_x(ii);
          uiy_1(ii)=beta*pos_y_11(ii)+alfa*spe_y_11(ii)+gamma*f_ct_y(ii)+omg*f_wall_y(ii);
      
          v_t(ii) = z_prey(ii,8,jj-1)*v_to*sigmf(z_prey(ii,5,jj-1),[10 0.5]);  % The magnitude of the velocity, a variable
%            v_t(ii)=v_to;  % the speed is a constant
          
          
          pix_1(ii)=v_t(ii)*uix_1(ii)/(sqrt(uix_1(ii)^2+uiy_1(ii)^2));
          piy_1(ii)=v_t(ii)*uiy_1(ii)/(sqrt(uix_1(ii)^2+uiy_1(ii)^2));
          
          % noise
          xita_prey = 2*noise_prey*pi*rand()-noise_prey*pi;
          pix_1(ii)=pix_1(ii)*cos(xita_prey)-piy_1(ii)*sin(xita_prey);
          piy_1(ii)=piy_1(ii)*cos(xita_prey)+pix_1(ii)*sin(xita_prey);
         
          gix_1(ii)=gix_1(ii)+pix_1(ii)*delta*prey_stop(ii);  % postion update, so each iteration should mention gix, giy and pix, piy
          giy_1(ii)=giy_1(ii)+piy_1(ii)*delta*prey_stop(ii);
   
    z_prey(ii,1,jj)=gix_1(ii);
    z_prey(ii,2,jj)=giy_1(ii);
    z_prey(ii,3,jj)=pix_1(ii);
    z_prey(ii,4,jj)=piy_1(ii);
    z_prey(ii,6:10,jj)=z_prey(ii,6:10,jj-1);  % keep unchanged unless any claim      
%% obtain theta_prey(ii)    
% dot_prey(ii) = z_prey(ii,3,jj-1)*pix_1(ii)+z_prey(ii,4,jj-1)*piy_1(ii);
% module_prey(ii) = sqrt((z_prey(ii,3,jj-1))^2+(z_prey(ii,4,jj-1))^2)*...
%     sqrt((pix_1(ii))^2+(piy_1(ii))^2);
% yuxuan_prey(ii)=dot_prey(ii)/module_prey(ii);
% theta_prey(ii) = acos(yuxuan_prey(ii));


dot_prey(ii)=pix_1(ii)-z_prey(ii,3,jj-1);
module_prey(ii)=sqrt((pix_1(ii)-z_prey(ii,3,jj-1))^2+(piy_1(ii)-z_prey(ii,4,jj-1))^2);
yuxuan_prey(ii)=dot_prey(ii)/module_prey(ii);            % cos(theta)=a/sart(a^2+b^2)
theta_prey(ii) = acos(yuxuan_prey(ii));



%% obtain fit, sbest and gbest/local best

if z_prey(ii,5,jj)==-1000 % obtain fit
    fit_prey(ii,jj)=0;
else 
    fit_prey(ii,jj)=fit_prey(ii,jj-1)+1;
end

if fit_prey(ii,jj)>sbest_prey(ii,6,jj-1) % obtain sbest
    sbest_prey(ii,6,jj)=fit_prey(ii,jj);
    sbest_prey(ii,1:5,jj)=z_prey(ii,6:10,jj);
else 
    sbest_prey(ii,:,jj)=sbest_prey(ii,:,jj-1);
    
end


for k = 1 :n_prey   %obtain local best
   d_11=sqrt((z_prey(k,1,jj-1)-z_prey(ii,1,jj-1))^2+...
              (z_prey(k,2,jj-1)-z_prey(ii,2,jj-1))^2);
          if d_11<R_nei_prey
              val = fit_prey(k,jj);
              fit_prey_neigh(ii,k) = val;
          else
              val = 0;
              fit_prey_neigh(ii,k) = val; 
          end

end

[lbest_prey(ii,6,jj),index]=max(fit_prey_neigh(ii,:));  % locl best

if lbest_prey(ii,6,jj)>lbest_prey(ii,6,jj-1)
    lbest_prey(ii,6,jj)=lbest_prey(ii,6,jj);
    lbest_prey(ii,1:5,jj)=z_prey(index,6:10,jj-1);
else
    lbest_prey(ii,6,jj)=lbest_prey(ii,6,jj-1);
    lbest_prey(ii,:,jj)=lbest_prey(ii,:,jj-1); 
end

%% energy update  and  rebirth of preys 

if z_prey(ii,5,jj)==-1000
    
 %-------------------------------------------------------------------   
%  % average tendency
%   AVG = (sum(z_prey(:,:,jj-1),1)-z_prey(ii,:,jj-1))/(n_prey-1);
%   z_prey(ii,6:10,jj) = AVG(6:10);   
 %------------------------------------------------------------------  
 % particle swarm optimization tendency
A = rand();B = rand(); C=A+B;

if rand()<0.1 % 10% Mutation rate
D=0.2*rand()-0.1;
else
    D=0;
end

   z_prey(ii,6,jj) = z_prey(ii,6,jj-1)+...
     cog*(A/C)*(sbest_prey(ii,1,jj)-z_prey(ii,6,jj-1))+...
     soc*(B/C)*(lbest_prey(ii,1,jj)-z_prey(ii,6,jj-1))+D;
   z_prey(ii,7,jj) = z_prey(ii,7,jj-1)+...
     cog*(A/C)*(sbest_prey(ii,2,jj)-z_prey(ii,7,jj-1))+...
     soc*(B/C)*(lbest_prey(ii,2,jj)-z_prey(ii,7,jj-1))+D;
   z_prey(ii,8,jj) = z_prey(ii,8,jj-1)+...
     cog*(A/C)*(sbest_prey(ii,3,jj)-z_prey(ii,8,jj-1))+...
     soc*(B/C)*(lbest_prey(ii,3,jj)-z_prey(ii,8,jj-1))+D;
   z_prey(ii,9,jj) = z_prey(ii,9,jj-1)+...
     cog*(A/C)*(sbest_prey(ii,4,jj)-z_prey(ii,9,jj-1))+...
     soc*(B/C)*(lbest_prey(ii,4,jj)-z_prey(ii,9,jj-1))+D;
   z_prey(ii,10,jj) = z_prey(ii,10,jj-1)+...
     cog*(A/C)*(sbest_prey(ii,5,jj)-z_prey(ii,10,jj-1))+...
     soc*(B/C)*(lbest_prey(ii,5,jj)-z_prey(ii,10,jj-1))+D;
 %-----------feasible area------------------------------------------
 if z_prey(ii,6,jj)<0
     z_prey(ii,6,jj)=0.1;
 else
 end
  if z_prey(ii,7,jj)<0
     z_prey(ii,7,jj)=0.1;
 else
  end
  if z_prey(ii,8,jj)<0
     z_prey(ii,8,jj)=0.1;
 else
  end
  if z_prey(ii,9,jj)<0
     z_prey(ii,9,jj)=0.1;
 else
  end
  if z_prey(ii,10,jj)<0
     z_prey(ii,10,jj)=0.1;
 else
 end
 %----------------------------------------------------------------
    z_prey(ii,5,jj)=0.8;   % the energy level of newbirth
    angle = 2*pi*rand();
    R_re = R_area*rand();
    a = R_re*cos(angle);
    b = R_re*sin(angle);
    c=zeros(1,2);
    z_prey(ii,1:4,jj)=[a,b,c];     
    gix_1(ii)=z_prey(ii,1,jj);
    giy_1(ii)=z_prey(ii,2,jj);
    pix_1(ii)=z_prey(ii,3,jj);
    piy_1(ii)=z_prey(ii,4,jj);
else 
    z_prey(ii,5,jj) = z_prey(ii,5,jj-1)+worse(ii)*k_escape*((v_t(ii)/v_to)^2+(theta_prey(ii)/pi)^2)+normal(ii)*U_c/Peer_num(ii); 
    z_prey(ii,6:10,jj)=z_prey(ii,6:10,jj-1); %update strategy sets
end  

%% normalization the energy of prey
if z_prey(ii,5,jj)<0 && z_prey(ii,5,jj)>-1000
    z_prey(ii,5,jj)=0;
elseif z_prey(ii,5,jj)>max_e
    z_prey(ii,5,jj)=max_e;
elseif z_prey(ii,5,jj)<min_e && z_prey(ii,5,jj)>0
    z_prey(ii,5,jj)=min_e;
end
    

    
             
   end
   
   %% data collection
   
   ESS_prey = mean(z_prey);
   ESS_pred = mean(z_pred);
   
   lif_prey = mean(fit_prey);
   max_prey = max(fit_prey);
   
     fprintf('times: %d \n',jj);
end

figure(1)

ESS_1=reshape(ESS_prey,10,timestep);


plot(ESS_1(5,:),'-m','LineWidth',2);
hold on
plot(ESS_1(6,:),'-g','LineWidth',2);
hold on
plot(ESS_1(7,:),'-k','LineWidth',2);
hold on
plot(ESS_1(8,:),'-b','LineWidth',2);
hold on
plot(ESS_1(9,:),'-r','LineWidth',2);
hold on
plot(ESS_1(10,:),'-y','LineWidth',2);
hold on
best_prey = ESS_1(:,end)
figure(2)

ESS_2=reshape(ESS_pred,10,timestep);


plot(ESS_2(5,:),'-m','LineWidth',2);
hold on
plot(ESS_2(6,:),'-g','LineWidth',2);
hold on
plot(ESS_2(7,:),'-k','LineWidth',2);
hold on
plot(ESS_2(8,:),'-b','LineWidth',2);
hold on
plot(ESS_2(9,:),'-r','LineWidth',2);
hold on
plot(ESS_2(10,:),'-y','LineWidth',2);
hold on


best_pred = ESS_2(:,end)

% figure (3)
%    plot(z_pred(:,1,1),z_pred(:,2,1),'^r');
%    hold on
%    plot(z_pred(:,1,end),z_pred(:,2,end),'^r','markerfacecolor','r');
%    hold on
%    
%    plot(z_prey(:,1,1),z_prey(:,2,1),'og');
%    hold on
%    plot(z_prey(:,1,end),z_prey(:,2,end),'og','markerfacecolor','g');
%    hold on
%    
%    axis equal
%  %=================================================      
%        
%         
%  aaa=z_pred(:,1,:);
%  bbb=z_pred(:,2,:);
%  aaa=reshape(aaa,n_pred,timestep);
%  bbb=reshape(bbb,n_pred,timestep);
%  for kkk=1:n_pred
%   plot(aaa(kkk,:),bbb(kkk,:),':b');
%   hold on
%  end
%  
%  ccc=z_prey(:,1,:);
%  ddd=z_prey(:,2,:);
%  ccc=reshape(ccc,n_prey,timestep);
%  ddd=reshape(ddd,n_prey,timestep);
%  for kkk=1:n_prey
%   plot(ccc(kkk,:),ddd(kkk,:),':k');
%   hold on
%  end   
 
 figure(4)
plot(lif_prey,'-m','LineWidth',2);
hold on  
plot(max_prey,'-g','LineWidth',2);
     