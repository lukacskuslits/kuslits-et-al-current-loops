%=============================Genetic algorithm=========================
% Solves the magnetic inverse problem in the outer core for m parameters
%-----------------------------------------------------------------------

function [eredmeny]=genalg_synth_plane_real_5(subname,fparam,subname_r,filename,looppar,dim,P,n,varargin)
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
p.addParamValue('norma','L2',@(x) ischar(x) && length(x)<30); %Set type of target function (norm) to call
p.addParamValue('childnr',2,@(x) isinteger(x) && x>0 && x<'n'); %Nr. of children
p.addParamValue('d',0.4,@(x) isscalar(x) && x>0 && x<1); %proportion of volumetric hypercube increase during crossover
%p.addParamValue('crnodes',2,@(x) isscalar(x) && x>0 && x<'m'); %Nr. of crossover points
p.addParamValue('pcr',0.6,@(x) isscalar(x) && x>0 && x<1); %Crossover probability
p.addParamValue('scale',1,@(x) isscalar(x) && x>0 && x<1000); %Parameter scaling
p.addParamValue('step',ones(140,1),@(x) inumeric(x)); %add resolution
p.addParamValue('FIX',141,@(x) isscalar(x) && x>0 && x<1000); %fix parameters from p nr fix to p nr looppar
p.addParamValue('brcr',300,@(x) isscalar (x) && x>0 && x<1000); %break-criterium
p.addParamValue('alfa',[],@(x) isnumeric(x)); % add alfa in case of Fletcher-Powell test



%Parse input arguments
p.parse(filename,looppar,dim,P,n,varargin{:});
p.Results

nloop=p.Results.nloop;
nii=p.Results.nii;
kiss=p.Results.kiss;
mutrate=p.Results.mutrate;
migrate=p.Results.migrate;
migtime=p.Results.migtime;
norma=p.Results.norma;
ngy=p.Results.childnr;
%crnodes=p.Results.crnodes;
d=p.Results.d;
pcr=p.Results.pcr;
scale=p.Results.scale;
step=p.Results.step;
FIX=p.Results.FIX;
brcr=p.Results.brcr;
alfa=p.Results.alfa;

filename=p.Results.filename;
%subname=p.Results.subname;
migtopo=p.Results.migtopo;



%  if length(varargin{:})>8
%     throw(exception)
%     disp('Tul sok bemeneti parameter!')
%  end
%=================================================================================================================
disp('Fajl betoltese.')


if dim>2
    load(filename,'-mat');
else
    ObjVal=load(filename,'ascii');
end
%subname=str2func(subname);
%subname=subname;
%--------------------------------------------------------------------------
%size(ObjVal)
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

if isempty(alfa)==1
xres=zeros(nloop,kiss);
yres=zeros(nloop,kiss);
zres=zeros(nloop,kiss);
Rres=zeros(nloop,kiss);
Ires=zeros(nloop,kiss);
tres=zeros(nloop,kiss);
lres=zeros(nloop,kiss);
end

pf=zeros(looppar,nloop,n,nii,P); % parameter space
p=zeros(m,n,nii,P); 
In=zeros(m,n,P); % initial parameter space
%size(p)

% setup parameter intervals
disp('Parameter tartomanyok szamitasa..')
partable=fparam(looppar,nloop);

pmax=zeros(m,1);
pmin=zeros(m,1);
% vec=zeros(dim+2,1);
for i=1:m
    pmax(i)=partable(i,1);
    pmin(i)=partable(i,2);
end

domx=pmax(1); % domain xsect width
domy=pmax(2);

%pmin

% for jj01=1:dim
%    vec(jj01)=pmax(jj01);    
% end
% 
% vec(dim+1)=n;
% vec(dim+2)=P;

% if isempty(alfa)==1
% fref=zeros(vec'); % measurement matrix
% else
% fref=zeros(1);
% end
%G=zeros(vec'); % model calculation matrix
%diff=zeros(vec'); % model-measurement difference matrix

eredmeny_control=zeros(looppar,nloop,nii/10+1,kiss);

for kis=1:kiss


%Initial parameter generation
%-------------------------------------------------------------------------
%pmin

disp('Elso nemzedek generalasa.')
%size(p)
 for nn=1:n
 for mm=1:m
if pmin(mm)<0 && pmin(mm)>-pi()
   pminn=pmin(mm);
   pmaxx=pmax(mm);
   stepp=step(mm);
 parfor pp=1:P 
 In(mm,nn,pp)=pminn+stepp*scale*randi([0 round((pmaxx-pminn)/(stepp*scale))]);
 p(mm,nn,1,pp)=In(mm,nn,pp);
 end
else
   pminn=pmin(mm);
   pmaxx=pmax(mm);
   stepp=step(mm);
 parfor pp=1:P
 In(mm,nn,pp)=pminn+stepp*randi([0 round((pmaxx-pminn)/(stepp))]);
 p(mm,nn,1,pp)=In(mm,nn,pp);
 end
end

 
 end
 end



for egyed0=1:n
for ii=1:nloop
for pp=1:P
       pf(:,ii,egyed0,1,pp)=In(looppar*ii-looppar+1:looppar*ii,egyed0,pp);  
end
end
end
% disp('Parameterek:')
% pf(:,:,:,1,4)


 if isempty(alfa)==1
 fref=log((ObjVal));% optimum function
 else
%   index=round(1/step*(pi+alfa/scale)'+1);
%   index=num2cell(index);
% 
%   R=struct;
%   R.type='()';
%   R.subs=index;
% 
%   fref=subsref(ObjVal,R);
    fref=min(min(min(min(min(min((ObjVal)))))));
 end

f=fref;
%f
%size(p)
%------------------------------------------------------------------------

target=zeros(n,nii,P);
fitness=zeros(n,nii,P);
SORTED=zeros(n,2,nii,P);

%Iteration
%================================================================================

for ni=1:nii-1 %iteration
 %size(p)
    %ni=1;

 
 %exclude stall or stick

   for popo=1:P
     if ni~=1 && min(target(:,ni,popo))-min(target(:,ni-1,popo))>=200
      p(:,:,ni,popo)=p(:,:,ni-1,popo);
      
     end
     if ni>5 && abs(min(target(:,ni,popo))-min(target(:,ni-4,popo)))<4....
     && abs(min(target(:,ni-1,popo))-min(target(:,ni-3,popo)))<3 && abs(min(target(:,ni-2,popo))-max(target(:,ni-1,popo)))<3
      p(:,:,ni,popo)=p(:,:,ni-5,popo);
      
     end
     
   end
 
   %Calculating forward problems and objective function
%-----------------------------------------------------------------------------------------------------------------
   disp('Direktfeladatok es a celfuggveny kiszamolasa.')
   for egyed1=1:n
   for ii=1:nloop
   for pp1=1:P
       pf(:,ii,egyed1,ni,pp1)=p(looppar*ii-looppar+1:looppar*ii,egyed1,ni,pp1);  
   end
   end
   end
   
 % exclude spatial homogeneity
if isempty(alfa)==1
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
end
 if ni>10  
  for sp=1:n
   for ps=1:n
   for popo=1:P
     pa=p(:,sp,ni,popo);
     pb=p(:,ps,ni-10,popo);

     if pa-pb==0
      pa=pmin; % constraint
      p(:,sp,ni,popo)=pa;
     end    
   end
   end
  end
 end
 
   for egyed1=1:n
   for ii=1:nloop
   parfor pp1=1:P
       pf(:,ii,egyed1,ni,pp1)=p(looppar*ii-looppar+1:looppar*ii,egyed1,ni,pp1);  
   end
   end
   end
 

% disp('Parameterek:')
% pf(:,:,:,ni,4)   

  try
   %size(p)
 
   for nn=1:n
   parfor popo=1:P
      
       
        target(nn,ni,popo)=subname(pf(:,:,nn,ni,popo),looppar,nloop,f,scale,domx,domy,norma,alfa);

       %G(:,:,nn,popo)=log(subname(pf(:,:,nn,ni,popo))); % ni'th generation
      
      
        %g=G(:,:,pop);
        %F=f(:,:,pop);
        
        %fitness(ni,pop)=1-max(max(abs(g-F)))/max(max(abs(F)));
   
        %fitness(ni,popsi)=OUT(ni,popsi); 
        
   end
   end
 
   %target

%%Fitness determination by linear ranking
%=========================================================================
   
   targ=zeros(n,2);
   
       for popsi=1:P
           for nn9=1:n
           targ(nn9,:)=[target(nn9,ni,popsi) nn9];
           end
           SORTED(:,:,ni,popsi)=sortrows(targ,-1);           
       end

    POST=zeros(n,P);   
   for egyed5=1:n
       for egyed6=1:n
       for pp9=1:P
       if SORTED(egyed5,2,ni,pp9)==egyed6
       POST(egyed6,pp9)=egyed5;
       end
       end
       end
   end

   
   ssp=3;
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
  
   %-------------------------------------------------------------------------------------
%Break-criteria
%----------------------------------------------------------------------------------------
  %size(p)
u1=0;  
    if mod(ni,10)==0 && ni>=10 % save control result in every 10th step
      MAX=min(min(target(:,ni,:),[],1),[],3);
      for egyed=1:n  
         %for lj=1:nii
             for pop3=1:P
               if target(egyed,ni,pop3)==MAX
                eredmeny_control(:,:,ni/10,kis)=pf(:,:,egyed,ni,pop3);
                save('eredmeny_control.mat','eredmeny_control');
                u1=1;
                break
               end
             end
              if u1==1;
              disp('Kontrol eredmeny mentese!')
              break 
%          end
%          %end
%         if u1==1
       
           % break
              end
      end           
    end
    
%disp('Parameterek:')
%pf(:,:,:,ni,4)
    %size(p)
    
    %Shutdown in case of reaching criterium
    u2=0;
    if min(min(target(:,ni,:),[],1),[],3)<=brcr
        MAX1=min(min(target,[],1),[],3);
       for egyed2=1:n 
        for pop2=1:P
            if target(egyed2,ni,pop2)==MAX1(1,ni)
             eredmeny(:,:,kis)=pf(:,:,egyed2,ni,pop2)
             target(egyed2,ni,pop2)
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
%c
%MAXIM
disp('Szelekcio.')
   % Roulette wheel selection
%------------------------------------------------------------------------------  
   psel=zeros(m,n,ni,P);
   sel=zeros(n,1);
   I=zeros(n,n,P);
   rws=zeros(P,1);
    fitc=zeros(n,P);
%constructing the wheel
 %size(p)
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
           for fixed=FIX:looppar
            if isempty(alfa)==1 && mod(ijj,fixed)~=0
            a=(0.5+0.5*randn(1,1)/3)*(-d)+(0.5+0.5*randn(1,1)/3)*(1+d);
            offspring(ijj,nno)=pany(ijj)*a+pap(ijj)*(1-a);
            else
            a=(0.5+0.5*randn(1,1)/3)*(-d)+(0.5+0.5*randn(1,1)/3)*(1+d);
            offspring(ijj,nno)=pany(ijj)*a+pap(ijj)*(1-a);
            end
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
      alfaloc=zeros(accur,1);
      alfaloc(randi(accur))=1;
      for mutvector=1:accur 
      delt(mutvector)=alfaloc(randi(accur))*2^(-mutvector);
      end
      if pmini(mut,popo2)<0
         sign=0;
      else
         sign=1;
      end
      Delta=sum(delt,1);
      pmut=p(mut,nnn,ni,popo2);
      eff=Mutmask(randi(m))*(pmaxi(mut,popo2)-sign*pmini(mut,popo2))*Mutshrink*Delta/2;
      
      for fixed2=FIX:looppar    
      if pmut+eff>=pmini(mut,popo2) && pmut+eff<=pmaxi(mut,popo2)&& mod(mut,fixed2)~=0
      pmut=pmut+eff;
      end
      p(mut,nnn,ni+1,popo2)=pmut;
      end
      end
  end
end
%size(p)
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
     %p(:,migi(migrant),ni+1,POPI)=zeros(m,1,1,1);
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
    % p(:,migi(migrant),ni+1,popo)=zeros(m,1,1,1);
     
     
  end   
 end
 migrant1=randi(n);
 migrant2=randi(n);
 p(:,migrant1,ni+1,1)=p(:,migrant2,ni+1,P);
 p(:,migrant2,ni+1,P)=zeros(m,1,1,1);
end
end

  
  MAXT=min(min(target,[],1),[],3);
     %draw fitness
     if mod(ni,1)==0
   
     hold on
     figure(kis), clf; 
     subplot(1,1,1);
     scatter(1:nii,MAXT,'ok')
     ylabel 'Elteres normaja'
     xlabel 'Iteracios lepes'
     title 'Aramhurkok keresese' 
     end
     drawnow;
disp('Lepesszam:')
ni
disp('Parameterek:')
pf(:,:,1,ni,3)
end
%==========================output structure==================================== 
disp('Az eredmenyek:')
eredmeny2=eredmeny(:,:,1)
save('eredmeny.mat','eredmeny2')
%save('eredmeny2.txt','eredmeny2','ascii')
%load('eredmeny.mat');

%% Calling linearised regularised inversion to produce final result
disp('Regularizalt javitas')
lam=0.2; %regularising factor
nit=4;  %linearised iteration step

res=regularised_synth3(eredmeny2,looppar,nloop,domx,domy,nit,lam,scale,ObjVal,subname_r);
eredmeny(:,:,1,kis)=res;
%=========================================================================

if isempty(alfa)==1
 xres(:,kis)=eredmeny(1,:,1,kis);
 yres(:,kis)=eredmeny(2,:,1,kis);
 zres(:,kis)=eredmeny(3,:,1,kis);
 Rres(:,kis)=eredmeny(4,:,1,kis);
 Ires(:,kis)=eredmeny(5,:,1,kis);
 tres(:,kis)=eredmeny(6,:,1,kis);
 lres(:,kis)=eredmeny(7,:,1,kis);
end
end   

%-------------------------------------------------------------------------------
if isempty(alfa)==1
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

filename_b = 'eredmeny_magnet.xlsx';

for kis=1:kiss
xlswrite(filename_b,res_b(:,:,kis),kis);
end
else
filename_b = 'eredmeny_fletcher_2.xlsx';
for kis=1:kiss
xlswrite(filename_b,eredmeny(:,1,kis),kis);
end
xlswrite(filename_b,alfa,kiss+1);
end
matlabpool close
toc
%*********************************************************************************
end


