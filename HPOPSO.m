clear all;
clc;
close all

format long;
global h Wssum Wcsum Vck Eck Time;

h=5;
WSSUM=[755640	1004640	962520	980040	1242840];
WCSUM=[512910	627940	632400	584160	841770];
ECK=[65.5	101.37	88.01	109.57	106.96];
Wssum=WSSUM(1:h);
Wcsum=WCSUM(1:h);
Eck=ECK(1:h);
Vck=100;
Time=[15	15	13	16	17];
Times=20;
function_name='F13';
[dim,CostFunction,ub, lb]  = Select_Functions(function_name);

%%%%%%%%PSO parameters
max_iteration=2000;         
swarm_size=30;           
particle_size=3;         
particle_min=[0 0.001*Wssum(1) 0.001*Wssum(1) ]; 
particle_max=[400 0.12*Wssum(1) 0.12*Wssum(1)];  
velocity_min=-1;            
velocity_max=1;           
c1=2;                     
c2=2;                    
wmax=1.0;                 
wmin=0.6;                  
CR=0.6;

BEST=zeros(Times,particle_size);
pingjun=zeros(max_iteration,particle_size+1);
for zong=1:Times
    zong
%%PSO
fitness_size=particle_size+1;                                      
particle=zeros(swarm_size,fitness_size,max_iteration);             
velocity=zeros(swarm_size,particle_size,max_iteration);           
pbest=zeros(swarm_size,fitness_size,max_iteration);               
gbest=zeros(1,fitness_size,max_iteration);                         
X1_new=zeros(swarm_size,particle_size);                            
value=zeros(1,swarm_size);                                         
 
z=1;
%%PSO initialization
for x=1:swarm_size
    for y=1:particle_size  
        particle(x,y,z)=(particle_min(y)+(particle_max(y)-particle_min(y))*rand);  
        velocity(x,y,(z+1))=rand;   
        %UpdateP
         particle(x,y,(z+1))=particle(x,y,z)+velocity(x,y,z+1);
             if particle(x,y,(z+1))>particle_max(y)
                  particle(x,y,(z+1))=particle_max(y)-(0.3e-10);
             end
             if particle(x,y,(z+1))<particle_min(y)
                  particle(x,y,(z+1))=particle_min(y)+(0.7e-10);
             end                 
    end
    particle(x,fitness_size,z)=fit(particle(x,:,z));  
end
pbest(:,:,z)=particle(:,:,z);
gbest(1,:,z)=pbest(1,:,z);  
for x=1:swarm_size
    if pbest(x,fitness_size,z)<gbest(1,fitness_size,z)
       gbest(1,:,z)=pbest(x,:,z);                            
    end
end


%%%%%%%%%%%%%%%%%%HPO initialization
  Convergence_curve = zeros(3,max_iteration);
B = 0.1;
 HPpos=rand(swarm_size,particle_size).*(ub-lb)+lb;
for i=1:size(HPpos,1)
HPposFitness(i)=CostFunction(HPpos(i,:));       
end
 [~,indx] = min(HPposFitness);
 Target = HPpos(indx,:);   % Target HPO
 TargetScore =HPposFitness(indx);
 Convergence_curve(1)=TargetScore;
 
for i=1:swarm_size
  X1_new(i,:,1)=particle(i,1:3,z);
end

%%%%%%%%%%%%%%%%%%%Iterative
for z=2:max_iteration
    %%%%%%%%%%
    for i=1:swarm_size
        if rand>CR   
            X1_new(i,:,z)=HPpos(i,:);
        else
            X1_new(i,:,z)=particle(i,1:3,z);
        end
    end
    
        %%%%%%%%%%%%%%%%%%%%%%select
 for i=1:swarm_size
     if fit(X1_new(i,:,z))<fit(X1_new(i,:,z-1))      
         X1_new(i,:,z)=X1_new(i,:,z);
     else 
        X1_new(i,:,z)=X1_new(i,:,z-1);
     end
 end
 for i=1:swarm_size 
     for j=1:particle_size
        particle(i,j,z)=X1_new(i,j,z);
        HPpos(i,j)=X1_new(i,j,z);
     end
 end
    
  %%%%%%%%%%%%%%HPO  
    
  c = 1 - z*((0.98)/max_iteration);   % Update C Parameter
    kbest=round(swarm_size*c);        % Update kbest
     for i = 1:swarm_size
            r1=rand(1,particle_size)<c;
            r2=rand;
            r3=rand(1,particle_size);
            idx=(r1==0);
            zz=r2.*idx+r3.*~idx;
        if rand<B
        xi=mean(HPpos);
        dist = pdist2(xi,HPpos);
        [~,idxsortdist]=sort(dist);
        SI=HPpos(idxsortdist(kbest),:);
        HPpos(i,:) =HPpos(i,:)+0.5*((2*(c)*zz.*SI-HPpos(i,:))+(2*(1-c)*zz.*xi-HPpos(i,:)));
        else
          for j=1:dim
            rr=-1+2*zz(j);
          HPpos(i,j)= 2*zz(j)*cos(2*pi*rr)*(Target(j)-HPpos(i,j))+Target(j);

          end
        end  
        HPpos(i,:) = min(max(HPpos(i,:),lb),ub);
        % Evaluation
        HPposFitness(i) = CostFunction(HPpos(i,:));
        % Update Target
        if HPposFitness(i)<TargetScore 
            Target = HPpos(i,:);
            TargetScore = HPposFitness(i);
        end
     end


 
% %Fitness  
    for x=1:swarm_size   
               particle(x,fitness_size,z)=fit(particle(x,:,z)); 
        if particle(x,fitness_size,z)<pbest(x,fitness_size,(z-1))
           pbest(x,:,z)=particle(x,:,z);
        else
           pbest(x,:,z)=pbest(x,:,(z-1));
        end   
    end
    
     gbest(1,:,z)=gbest(1,:,(z-1));
        for x=1:swarm_size
            if pbest(x,fitness_size,z)<gbest(1,fitness_size,(z-1))
               gbest(1,:,z)=pbest(x,:,z); 
            end
        end
    
%Update
    for x=1:swarm_size
        for y=1:particle_size
          %UpdateV
           velocity(x,y,(z+1))=(wmax-(z/max_iteration)*(wmax-wmin))*velocity(x,y,z)+c1*rand*(pbest(x,y,z)-particle(x,y,z))+c2*rand*(gbest(1,y,z)-particle(x,y,z));
              if velocity(x,y,(z+1))>velocity_max
                    velocity(x,y,(z+1))=velocity_max-(0.3e-10);
              end
              if velocity(x,y,z+1)<velocity_min
                    velocity(x,y,z+1)=velocity_min+(0.7e-10);
              end
            %UpdateP
             particle(x,y,(z+1))=particle(x,y,z)+velocity(x,y,z+1);
             if particle(x,y,(z+1))>particle_max(y)
                  particle(x,y,(z+1))=particle_max(y)-(0.3e-10);
             end
             if particle(x,y,(z+1))<particle_min(y)
                  particle(x,y,(z+1))=particle_min(y)+(0.7e-10);
             end               
        end
    end
end


%%%%%Results
BEST(zong,1:4)=gbest(1,1:4,max_iteration);
%%%%%Fitness function curve
for i=1:max_iteration
    for j=1:4
        pingjun(i,j)= pingjun(i,j)+gbest(1,j,i)/Times;
    end
end
end
AVE=mean(BEST(:,4));