%=============================Genetic algorithm=========================
% Solves the magnetic inverse problem in the outer core for m parameters
%-----------------------------------------------------------------------

function [eredmeny]=genalg_synth_plane_real_2(subname,fparam,filename,looppar,dim,P,n,varargin)
matlabpool open;

%%**********************************************************************
%                               Instructions
%
% use @funct_name in command to define parameter space and forward function
% to operate..
%
% essential arguments for funct. genalg are: subnamem, fparam (see above),
% values looppar, dim, P, and n..


clc
clf
tic

disp 'Vajon hol milyen aramgyuruk lehetnek a magban es mennyi??'

%Parse input, algorithm variables
%-----------------------------------------------------------------------------------------------------------------
p=inputParser;
p.addRequired('filename',@(x) ischar (x) && length(x)<30); %name of objective funct. file to load
%p.addRequired('subname',@(x) ischar (x) && length(x)<30); %name of forward funct. subroutine to call
%p.addRequired('m',@(x) isscalar(x) && x>1 && x<10000); % nr of params
p.addRequired('looppar',@(x) isscalar(x) && x>0 && x<10); %Nr. of features per source (in case of source identification probl.) 
%..or nr. of parameters to estimate
p.addRequired('dim',@(x) isscalar(x) && x>0 && x<10); % dimension of the phenotype function of the problem
p.addRequired('P',@(x) isscalar(x) && x>1 && x<200); % nr of populations
p.addRequired('n',@(x) isscalar(x) && x>1 && x<1000); % nr of specimen 
p.addParamValue('nloop',1,@(x) isscalar(x) && x>0 && x<100); %Nr. of sources (in case of source identification probl.)
p.addParamValue('nii',10,@(x) isscalar(x) && x>1 && x<500); %Iteration number
p.addParamValue('kiss',1,@(x) isscalar(x) && x>0 && x<5); %Trial number
p.addParamValue('mutrate',0.1,@(x) isscalar(x) && x>0 && x<1); % Probability (rate) of mutation
p.addParamValue('migrate',0.1,@(x) isscalar(x) && x>0 && x<1); % Probability (rate) of migration between populations
p.addParamValue('migtopo',@(x) ischar(x) && length(x)<30); %Set migration topology (ring or nonrestricted)
p.addParamValue('migtime',@(x) isinteger(x) && x>0 && x<30); %Set migration time window
p.addParamValue('childnr',2,@(x) isinteger(x) && x>0 && x<'n'); %Nr. of children
p.addParamValue('d',0.4,@(x) isscalar(x) && x>0 && x<1); %proportion of volumetric hypercube increase during crossover
%p.addParamValue('crnodes',2,@(x) isscalar(x) && x>0 && x<'m'); %Nr. of crossover points
p.addParamValue('pcr',0.6,@(x) isscalar(x) && x>0 && x<1); %Crossover probability


%Parse input arguments
p.parse(filename,looppar,dim,P,n,varargin{:});
p.Results

nii=p.Results.nii;
kiss=p.Results.kiss;
mutrate=p.Results.mutrate;
migrate=p.Results.migrate;
ngy=p.Results.childnr;
%crnodes=p.Results.crnodes;
pcr=p.Results.pcr;
looppar=p.Results.looppar;
nloop=p.Results.nloop;
migtime=p.Results.migtime;
d=p.Results.d;
filename=p.Results.filename;
%subname=p.Results.subname;
migtopo=p.Results.migtopo;

%  if length(varargin{:})>8
%     throw(exception)
%     disp('Tul sok bemeneti parameter!')
%  end
%=================================================================================================================
disp('Fajl betoltese.')
fielddata_plane_synth_1=load(filename,'ascii');
%subname=str2func(subname);
%subname=subname;
%-----------------------------------------------------------------------

% P=6;
% n=16;
% k=41;
% Zend=25;
% nii=80;
% kis=1;
% looppar=7;
% nloop=20;
m=looppar*nloop;
% pcr=0.6;
% load('fielddata_plane_synth_2.txt','ascii');
% ngy=2;
% mutrate=15/m;
% migrate=0.05;

%Result structure array
eredmeny=zeros(looppar,nloop,kiss);
xres=zeros(nloop,kiss);
yres=zeros(nloop,kiss);
zres=zeros(nloop,kiss);
Rres=zeros(nloop,kiss);
Ires=zeros(nloop,kiss);
tres=zeros(nloop,kiss);
lres=zeros(nloop,kiss);


pf=zeros(looppar,nloop,n,nii,P); % parameter space
p=zeros(m,n,nii,P); 
In=zeros(m,n,P); % initial parameter space
size(p)

% setup parameter intervals
disp('Parameter tartomanyok szamitasa..')
partable=fparam(looppar,nloop);

pmax=zeros(m,1);
pmin=zeros(m,1);
vec=zeros(dim+2,1);
for i=1:m
    pmax(i)=partable(i,1);
    pmin(i)=partable(i,2);
end

for jj01=1:dim
   vec(jj01)=pmax(jj01);    
end

vec(dim+1)=n;
vec(dim+2)=P;

f=ones(vec'); % measurement matrix
G=zeros(vec'); % model calculation matrix
diff=zeros(vec'); % model-measurement difference matrix



for kis=1:kiss


%Initial parameter generation
%------------------------------------------------------------------------

disp('Elso nemzedek generalasa.')
size(p)
for nn=1:n
for mm=1:m
parfor pp=1:P
In(mm,nn,pp)=pmin(mm)+(0.5+0.5*randn(1,1)/3)*(pmax(mm)-pmin(mm));
p(mm,nn,1,pp)=In(mm,nn,pp);
end
end
end


for egyed0=1:n
for ii=1:nloop
parfor pp=1:P
       pf(:,ii,egyed0,1,pp)=In(looppar*ii-looppar+1:looppar*ii,egyed0,pp);  
end
end
end
disp('Parameterek:')
pf(:,:,1,1,4)

for nnn=1:n
parfor popo=1:P
 
f(:,:,nnn,popo)=log((fielddata_plane_synth_1));% optimum function

end
end
size(p)
%------------------------------------------------------------------------

target=zeros(n,nii,P);
fitness=zeros(n,nii,P);


%Iteration
%================================================================================

for ni=1:nii-1 %iteration
 size(p)
    %ni=1;

 
 %exclude stall or stick

   for popo=1:P
%      if ni~=1 && max(fitness(:,ni,popo))-max(fitness(:,ni-1,popo))<=-1
%       p(:,:,ni,popo)=p(:,:,ni-1,popo);
%       
%      end
     if ni>5 && abs(max(target(:,ni,popo))-max(target(:,ni-4,popo)))<0.01...
     && abs(max(target(:,ni-1,popo))-max(target(:,ni-3,popo)))<1e-5 && abs(max(target(:,ni-2,popo))-max(target(:,ni-1,popo)))<1e-5
      p(:,:,ni,popo)=p(:,:,ni-5,popo);
      
     end 
   end
 
   %Calculating forward problems and objective function
%-----------------------------------------------------------------------------------------------------------------
   disp('Direktfeladatok es a celfuggveny kiszamolasa.')
   for egyed1=1:n
   for ii=1:nloop
   parfor pp1=1:P
       pf(:,ii,egyed1,ni,pp1)=p(looppar*ii-looppar+1:looppar*ii,egyed1,ni,pp1);  
   end
   end
   end
 % exclude spatial homogeneity
  for egyed7=1:n
   for popo=1:P
     p3=pf(1:3,:,egyed7,ni,popo);
     p5=pf(1:5,:,egyed7,ni,popo);
     if p3==0
         egyed7, popo
     end
     
     for sp=1:nloop
       for ps=nloop:1
        if ps~=sp && p3(sp)-p3(ps)<0.5
            p5(ps)=0; % constraint
            %p(:,:,ni-2,pop);
        end
       end
     end
     pf(1:5,:,egyed7,ni,popo)=p5;
   end
  end

disp('Parameterek:')
pf(:,:,1,ni,4)   

  try
   size(p)
   for nn=1:n
   parfor popsi=1:P

      G(:,:,nn,popsi)=log(subname(pf(:,:,nn,ni,popsi))); % ni'th generation
      
      
        %g=G(:,:,pop);
        %F=f(:,:,pop);
        
        %fitness(ni,pop)=1-max(max(abs(g-F)))/max(max(abs(F)));
   
        %fitness(ni,popsi)=OUT(ni,popsi); 
        diff(:,:,nn,popsi)=G(:,:,nn,popsi)-f(:,:,nn,popsi);
       target(nn,ni,popsi)=1/abs(norm(diff(:,:,nn,popsi)));
   end
   end
   

%%Fitness determination by linear ranking
%=========================================================================
   SORT=zeros(n,2,nii,P);
   targ=zeros(n,2);
   
       for popsi=1:P
           for nn9=1:n
           targ(nn9,:)=[target(nn9,ni,popsi) nn9];
           end
           SORT(:,:,ni,popsi)=sortrows(targ,1);           
       end
     
       POST=zeros(n,P);
   for egyed5=1:n
       for egyed6=1:n
       for pp9=1:P
       if SORT(egyed5,2,ni,pp9)==egyed6
       POST(egyed6,pp9)=egyed5;
       end
       end
       end
   end
   
   ssp=2;
   for nn4=1:n
       parfor popsi=1:P
        fitness(nn4,ni,popsi)=2-ssp+2*(ssp-1)*((POST(nn4,popsi))-1)/(n-1);
       end
   end   
   
%=========================================================================
%        parfor popsi1=1:P
%            SUMT(ni,popsi1)=sum(target(:,ni,popsi1),1);           
%        end
%   
% 
%    for nn4=1:n
%        parfor popsi=1:P
%         fitness(nn4,ni,popsi)=target(nn4,ni,popsi)/SUMT(ni,popsi);
%        end
%    end

     for nn5=1:n
       parfor popsi2=1:P   

        if isnan(fitness(nn5,ni,popsi2))==1
            fitness(nn5,ni,popsi2)=0.01;
            %disp('nan!')
        end
        if fitness(nn5,ni,popsi2)<=0
         fitness(nn5,ni,popsi2)=0.01;
         %disp('neG!')
        end
       end
     end

   
   catch Me
       fclose('all');
       disp(Me)
       disp('Hiba tortent a celfuggveny szamitasakor!')
  end
  G(:,:,1,1)
   %-------------------------------------------------------------------------------------
%Break-criteria
%----------------------------------------------------------------------------------------
  size(p)
u1=0;  
    if max(max(target(:,ni,:),[],1),[],3)<.1 && ni>=nii-1 % shutdown iteration, and end current trial in case of permanent stick
      MAX=max(max(max(target)));
      for egyed=1:n  
         for lj=1:ni
             for pop3=1:P
             if fitness(egyed,lj,pop3)==MAX
             eredmeny(:,:,kis)=pf(:,:,egyed,lj,pop3);
             disp('Sajnos beragadt!')
             u1=1;
             break
             end
             end
         if u1==1;
             break 
         end
         end
        if u1==1
            break
        end
      end
      break     
    end
    
disp('Parameterek:')
pf(:,:,1,ni,4)
    size(p)
    
    %Shutdown in case of reaching criterium
    u2=0;
    if max(max(target(:,ni,:),[],1),[],3)==.1
        MAX1=max(max(target,[],1),[],3);
       for egyed2=1:n 
        for pop2=1:P
            if target(egyed2,ni,pop2)==MAX1(1,ni)
             eredmeny(:,:,kis)=pf(:,:,egyed2,ni,pop2);
             u2=1;
             break
            end
        end
        if u2==1
            break 
        end
       end
        disp('Elerte a konvergencia-kriteriumot:')
        disp('ni')
        break   
    end
disp('Parameterek:')
pf(:,:,1,ni,4)
    
 %----------------------------------------------------------------------------------------------------
     MAXIM=max(max(fitness(:,ni,:),[],1),[],3);
         if 1<MAXIM && MAXIM<=10
            c=100;
         end
         if 0.1<MAXIM && MAXIM<=1
            c=1000;
         end
         if MAXIM<=0.1
            c=10000;
         end 

disp('Szelekcio.')
   % Roulette wheel selection
%------------------------------------------------------------------------------  
   psel=zeros(m,n,ni,P);
   sel=zeros(n,1);
   I=zeros(n,n,P);
   rws=zeros(P,1);
    fitc=zeros(n,P);
%constructing the wheel
 size(p)
   for pop5=1:P
    fitc(:,pop5)=c*(fitness(:,ni,pop5));
    if fitc(:,pop5)==0
        fitc(:,pop5)=0.01;
    end
    rws(pop5)=sum(round(fitc(:,pop5)),1);
    RW=zeros(rws(pop5),1);
    
    for l=1:n
    I(1:l,l,pop5)=round(fitc(1:l,pop5));
    
    j=sum(I,1);
    
       
    if l==1
        
     for  kw=1:j(1,1,pop5); 
         RW(kw)=1;
     end 
      
    else
     for  kw=j(1,l-1,pop5)+1:j(1,l,pop5);  
        RW(kw)=l;
     end
    end


    end
    for e=1:n
        sel(e)=randi(rws(pop5));
        psel(:,e,ni,pop5)=p(:,RW(sel(e)),ni,pop5); 
    end
   end


    p(:,:,ni,:)=psel(:,:,ni,:);
%---------------------------------------------------------------------------------------------

 %crossover and mutation
%--------------------------------------------------------------------------------------------- 

 % initializing crossover and mutation effect
     nmut=round(m*mutrate); 
     pmaxi=zeros(m,P);
     pmini=zeros(m,P);
     sign=[1;-1];
     accur=20;
     Mutshrink=1;
      % probabilistic determination of mutation effect 
     crossover=[ones(round(n*pcr),1);zeros(n-round(n*pcr),1)]; % crossover probability is set to pcr%
     % mutation probability is set to mutrate
     Mutmask=zeros(m,1);
     for maskelement=1:nmut
         Mutmask(maskelement)=sign(randi(2)); 
     end

     
            for segg=1:P
            pmaxi(:,segg)=round(pmax(:));
            pmini(:,segg)=round(pmin(:));
            end           
            

 %-----------------------------------------------------------------------
 disp('Rekombinacio.')
for popo=1:P
 

% intermediate recombination
%%-------------------------------------------------------------------
offspring=zeros(m,ngy);
for spsp=1:n/2
any=randi(n);
ap=randi(n-1);
pany=p(:,any,ni,popo);

papsel=cat(2,p(:,1:any-1,ni,popo),p(:,any+1:n,ni,popo));
pap=papsel(:,ap);

if crossover(randi(n))==1
for nno=1:ngy
    for ijj=1:m
    offspring(ijj,nno)=pany(ijj);
    if mod(ijj,6)~=0 && mod(ijj,7)~=0
    a=(0.5+0.5*randn(1,1)/3)*(-d)+(0.5+0.5*randn(1,1)/3)*(1+d);
    offspring(ijj,nno)=pany(ijj)*a+pap(ijj)*(1-a);
    end
    end
  
    p(:,2*spsp+nno-2,ni+1,popo)=offspring(:,nno);
end
else 
p(:,any,ni+1,popo)=pany;
    if ap<any
       p(:,ap,ni+1,popo)=pap;
    else    
       p(:,ap+1,ni+1,popo)=pap;
    end
end
end
end
%%--------------------------------------------------------------------

 %----------------------------------------------------------------------------------
  disp('Mutacio.')

  delt=zeros(accur,1);
  
for popo2=1:P 
  for nnn=1:n
      for mut=1:m
      alfa=zeros(accur,1);
      alfa(randi(accur))=1;
      for mutvector=1:accur 
      delt(mutvector)=alfa(randi(accur))*2^(-mutvector);
      end
      Delta=sum(delt,1);
      pmut=p(mut,nnn,ni,popo2);
      eff=Mutmask(randi(m))*(pmaxi(mut,popo2)-pmini(mut,popo2))*Mutshrink*Delta/2;
      
      if pmut+eff>=pmini(mut,popo2) && pmut+eff<=pmaxi(mut,popo2)
      pmut=pmut+eff;
      end
      p(mut,nnn,ni+1,popo2)=pmut;
      end
  end
end
size(p)
%% Reinsertion
%%====================================================================
% Random reinsertion
disp('Reinzercio')

shuffle=randperm(n);

for nn8=1:n
    for popo3=1:P
    p(:,nn8,ni+1,popo3)=p(:,shuffle(1,nn8),ni+1,popo3);
    end
end
  
%--------------------------------------------------------------------------------------------------------------------

%migration
%--------------------------------------------------------------------------------------------------------------------
%uniformly distributed migrate% migration between the populations


migrants=round(n*migrate);
mige=zeros(migrants,1);
migi=zeros(migrants,1);

%migration in nonrestricted topology
if strcmp(migtopo,'rand')==1
 for popo=1:P
     POPI=randi(P);
     POPE=randi(P);

     for migrant=1:migrants
      mige(migrant)=randi(n);
      migi(migrant)=randi(n);
      
     p(:,mige(migrant),ni+1,POPE)=p(:,migi(migrant),ni+1,POPI);
     p(:,migi(migrant),ni+1,POPI)=zeros(m,1,1,1);
     end  
     
 end
end

%migration in ring topology
if strcmp(migtopo,'ring')==1
if mod(ni,migtime)==0
    disp('Migracio.')
 for migrant=1:migrants
  for popo=1:P-1
     
      mige(migrant)=randi(n);
      migi(migrant)=randi(n);
     
     p(:,mige(migrant),ni+1,popo+1)=p(:,migi(migrant),ni+1,popo);
     p(:,migi(migrant),ni+1,popo)=zeros(m,1,1,1);
     
     
  end   
 end
 migrant1=randi(n);
 migrant2=randi(n);
 p(:,migrant1,ni+1,1)=p(:,migrant2,ni+1,P);
 p(:,migrant2,ni+1,P)=zeros(m,1,1,1);
end
end

  
  MAXT=max(max(target,[],1),[],3);
     %draw fitness
     if mod(ni,1)==0
   
     hold on
     figure(kis), clf; 
     subplot(1,1,1);
     scatter(1:nii,MAXT,'ok')
     ylabel 'Korrelacio'
     xlabel 'Iteracios lepes'
     title 'Aramhurkok keresese' 
     end
     drawnow;
disp('Lepesszam:')
ni
disp('Parameterek:')
pf(:,:,1,ni,4)
end
%==========================output structure==================================== 
disp('Az eredmenyek:')


 
 xres(:,kis)=eredmeny(1,:,1,kis);
 yres(:,kis)=eredmeny(2,:,1,kis);
 zres(:,kis)=eredmeny(3,:,1,kis);
 Rres(:,kis)=eredmeny(4,:,1,kis);
 Ires(:,kis)=eredmeny(5,:,1,kis);
 tres(:,kis)=eredmeny(6,:,1,kis);
 lres(:,kis)=eredmeny(7,:,1,kis);

end   

%-------------------------------------------------------------------------------
res=struct('x',{xres},'y',{yres}...
 ,'z',{zres},'R',{Rres},'I',{Ires},'t',{tres},...
 'l',{lres});

disp('x:')
res.xres
disp('y:')
res.yres
disp('z:')
res.zres
disp('R:')
res.Rres
disp('I:')
res.Ires
disp('t:')
res.tres
disp('l:')
res.lres

res_b(1,:,:)=res.xres;
res_b(2,:,:)=res.yres;
res_b(3,:,:)=res.zres;
res_b(4,:,:)=res.Rres;
res_b(5,:,:)=res.Ires;
res_b(6,:,:)=res.tres;
res_b(7,:,:)=res.lres;

filename_b = 'eredmeny_c.xlsx';

for kis=1:kiss
xlswrite(filename_b,res_b(:,:,kis),kis);
end
toc
%*********************************************************************************
end


