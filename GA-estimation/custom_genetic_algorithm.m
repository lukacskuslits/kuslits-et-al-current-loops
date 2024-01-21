%=============================Genetic algorithm=========================
% Solves the magnetic inverse problem in the outer core for m parameters
%-----------------------------------------------------------------------

function [eredmeny]=custom_genetic_algorithm(subname,fparam,table,filename,looppar,dim,P,n,domx,domy, coef, pow, varargin)

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
p.addRequired('looppar',@(x) isscalar(x) && x>0 && x<20); %Nr. of features per source (in case of source identification probl.) 
%..or nr. of parameters to estimate
p.addRequired('dim',@(x) isscalar(x) && x>0 && x<10); % dimension of the phenotype function of the problem
p.addRequired('P',@(x) isscalar(x) && x>1 && x<200); % nr of populations
p.addRequired('n',@(x) isscalar(x) && x>1 && x<1000); % nr of specimen 
p.addParamValue('nloop',1,@(x) isscalar(x) && x>0 && x<500); %Nr. of sources (in case of source identification probl.)
p.addParamValue('nii',10,@(x) isscalar(x) && x>1 && x<3000); %Iteration number
p.addParamValue('kiss',1,@(x) isscalar(x) && x>0 && x<5); %Trial number
p.addParamValue('mutrate',0.1,@(x) isscalar(x) && x>0 && x<1); % Probability (rate) of mutation
p.addParamValue('migrate',0.1,@(x) isscalar(x) && x>0 && x<1); % Probability (rate) of migration between populations
p.addParamValue('migtopo','rand',@(x) ischar(x) && length(x)<30); %Set migration topology (ring or nonrestricted)
p.addParamValue('migtime',5,@(x) isscalar(x) && x>0 && x<30); %Set migration time window
p.addParamValue('norma','L2',@(x) ischar(x) && length(x)<30); %Set type of target function (norm) to call
p.addParamValue('childnr',2,@(x) isinteger(x) && x>0 && x<'n'); %Nr. of children
p.addParamValue('d',0.4,@(x) isscalar(x) && x>0 && x<1); %proportion of volumetric hypercube increase during crossover
%p.addParamValue('crnodes',2,@(x) isscalar(x) && x>0 && x<'m'); %Nr. of crossover points
p.addParamValue('pcr',0.6,@(x) isscalar(x) && x>0 && x<1); %Crossover probability
p.addParamValue('scale',1,@(x) isscalar(x) && x>0 && x<1000); %Parameter scaling
p.addParamValue('const',1,@(x) isscalar(x) && x>-1e6 && x<1e6); %Parameter multiplication
p.addParamValue('step',ones(2000,1),@(x) inumeric(x)); %add resolution
p.addParamValue('free','',@(x) ischar(x) && length(x)<30); %fix parameters from p nr fix to p nr looppar
p.addParamValue('free2',[],@(x) isnumeric(x)); %fix parameters of each source between 1 & nloop
p.addParamValue('flmode','f',@(x) ischar(x) && length(x)<30); %fix parameters from p nr fix to p nr looppar
%p.addParamValue('sPOP',10,@(x) isscalar(x) && x>0 && x<50); %Nr. of stored results
p.addParamValue('brcr',1e-9, @(x) isscalar(x) && x>0 && x<1000); %break criterium
p.addParamValue('alfa',[],@(x) isnumeric(x)); % add alfa in case of Fletcher-Powell test
p.addParamValue('saveAs',[],@(x) ischar(x)); % add alfa in case of Fletcher-Powell test
p.addParamValue('lconfig',@(x) isscalar(x) && x>1 && x<41); % lounch configuration
p.addParamValue('broken', 0,@(x) isscalar(x) && x>=0 && x<=1); % start from previous step if crashed


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
flmode=p.Results.flmode;
%sPOP=p.Results.sPOP;
brcr=p.Results.brcr;
alfa=p.Results.alfa;
saveAs=p.Results.saveAs;
lconfig=p.Results.lconfig;
broken = p.Results.broken;

filename=p.Results.filename;
migtopo=p.Results.migtopo;
%matlabpool open

parpool(lconfig);

save_interval = 1;
%=================================================================================================================
disp('Fajl betoltese.')


if dim>2
    ObjVal=load(filename,'-mat');
    ObjVal=ObjVal.f;
else
    ObjVal=load(filename,'-ascii');
end

if length(ObjVal(:,1))==1
ObjVal=ObjVal';
end
%--------------------------------------------------------------------------
m=looppar*nloop;


%Result structure array
%eredmeny=zeros(looppar,nloop,kiss);

if isempty(alfa)==1
xres=zeros(nloop,kiss);
yres=zeros(nloop,kiss);
zres=zeros(nloop,kiss);
Rres=zeros(nloop,kiss);
Ires=zeros(nloop,kiss);
tres=zeros(nloop,kiss);
lres=zeros(nloop,kiss);
end


pf=zeros(looppar,nloop,n); % parameter space
p0=zeros(m,n,P); 
In=zeros(m,n,P); % initial parameter space
%size(p)
% if free=='e'
%     free=['t','l','d','a','R','I','m','z','h'];
% end



%Array of flags assigned to data types
   flags=zeros(looppar,nloop);
   
   if flmode=='i'
       for fl=1:looppar
           flags(fl,:)='i';
       end
   elseif flmode=='f'
       for fl=1:looppar
           flags(fl,:)='f';
       end
   else
       for fl=1:nloop
           flags(:,fl)=flmode';
       end
   end
   size(flags)
   for jk=1:nloop
   flag(looppar*jk-looppar+1:looppar*jk)=flags(:,jk);
   end
   flag=flag';
% setup parameter intervals
disp('Parameter tartomanyok szamitasa..')
%table
partable=fparam(table);

pmax=zeros(m,1);
pmin=zeros(m,1);

for i=1:m
    pmax(i)=partable(i,1);
    pmin(i)=partable(i,2);
end
disp(pmax),disp(pmin)
%if borken == 0
eredmeny_control=zeros(looppar,nloop,nii/save_interval+1,kiss);
%else

%load('eredmeny_control_oc_h2.mat'); 
%end
%POPNI=zeros(sPOP,nii);
PlotPop=zeros(nii,1);
for kis=1:kiss


%Initial parameter generation
%-------------------------------------------------------------------------

if broken == 0
disp('Elso nemzedek generalasa.')
if max(pmax)>15
 for nn=1:n
 for mm=1:m
if pmin(mm)<0 && pmin(mm)>-pi()
   pminn=pmin(mm);
   pmaxx=pmax(mm);
   stepp=step(mm);
 for pp=1:P 
 In(mm,nn,pp)=pminn+stepp*scale*randi([0 round((pmaxx-pminn)/(stepp*scale))]);
 p0(mm,nn,pp)=In(mm,nn,pp);
 end
else

   sign=[-1 1];
   if mod(mm,6)==4
   sign=sign(randi(2));
   else
   sign=1;
   end
   pminn=pmin(mm);
   pmaxx=pmax(mm);

   stepp=step(mm);
 for pp=1:P
 In(mm,nn,pp)=sign*pminn+stepp*sign*randi([0 round((pmaxx-pminn)/(stepp))]);
 p0(mm,nn,pp)=In(mm,nn,pp);
 end
end

 
 end
 end
else
 for nn=1:n
 for mm=1:m  
   sign=1;
   pminn=pmin(mm);
   pmaxx=pmax(mm);

   stepp=step(mm);
 for pp=1:P
 In(mm,nn,pp)=sign*pminn+.1*stepp*sign*randi([0 round(10*(pmaxx-pminn)/(stepp))]);
 p0(mm,nn,pp)=In(mm,nn,pp);
 end
 end
 end
end
else
    disp('Folytat!')
    load('P_loop_total.mat')
    p0 = P_1;
end
 
% disp('Parameterek:')
% pf(:,:,:,1,4)


 if isempty(alfa)==1
 fref=ObjVal;% optimum function
 else
    fref=min(min(min(min(min(min((ObjVal)))))));
 end
 
f=fref;
%[f,indices,bounds]=drop_ind(f);
%w=1;%weights_ready(f,domx,domy);

%argin=struct('p_ini',p_ini,'domx',domx,'domy',domy,'const',const,'free1',free,'free2',free2,'indices',indices,'bounds',bounds);
%argin=struct('domx',domx,'domy',domy);

%f
%size(p)
%------------------------------------------------------------------------
plast5=zeros(m,n,P);
pdisp=zeros(looppar,nloop,n,P);
target=zeros(n,nii,P);
fitness=zeros(n,1);
DT=zeros(nii,1);
MUT=zeros(nii,1);
thresholdD=0.05;
thresholdM=7;
%Iteration
%================================================================================

for ni=1:nii-1 %iteration
        d1=d;
        mutrate1=mutrate;
%Alternatively        
%
%     if ni==1
%         d1=d;
%         mutrate1=mutrate;
%     else
%      d1=d*abs((log(100*brcr+1)-abs(log(100*MAXT(ni)+1)))/(log(100*brcr+1)));
%      mutrate1=mutrate*abs((log(100*brcr+1)-abs(log(100*MAXT(ni)+1)))/(log(100*brcr+1)));
%     end
    
for popo=1:P 
    
    p=p0(:,:,popo);
 %p(:,:,ni,2)
 
 %exclude stall or stick

   
     if ni>5 && abs(min(target(:,ni,popo))-min(target(:,ni-4,popo)))<1e-15...
     && abs(min(target(:,ni-1,popo))-min(target(:,ni-3,popo)))<1e-15 && abs(min(target(:,ni-2,popo))-max(target(:,ni-1,popo)))<1e-15
      p=plast5(:,:,popo);
      
     end
     
   
 
   %Calculating forward problems and objective function
%-----------------------------------------------------------------------------------------------------------------
   if popo==1
    disp('Direktfeladatok es a celfuggveny kiszamolasa.')
   end
   for egyed1=1:n
   for ii=1:nloop
       pf(:,ii,egyed1)=p(looppar*ii-looppar+1:looppar*ii,egyed1);  
   end
       pdisp(:,:,egyed1,popo)=pf(:,:,egyed1);
   end
try
   if dim>2
       f1 = squeeze(f(1,:,:));
       f2 = squeeze(f(2,:,:));
   end
    
   parfor nn=1:n
        pcal=reshape(p(:,nn),looppar,nloop);
        if dim>2
            [res_field,res_sv]=subname(pcal,domx,domy,coef,pow);
        else
            res=subname(pcal,domx,domy,coef,pow);
        end
        if strcmp(norma,'rel')==1
            if dim<=2
               target(nn,ni,popo)=sum(sum((abs(res-f))))/(domx*domy);
            else
               target(nn,ni,popo)=sum(sum((abs(res_field-f1))))/(domx*domy)+sum(sum((abs(res_sv-f2))))/(10*domx*domy);
            end                  
        end
   end

disp('Celfuggveny')
targni = target(:,ni,:);
best_guess_ni = min(min(targni(targni>0)));
disp(best_guess_ni)
disp(targni(targni==best_guess_ni))
   %target

%%Fitness determination by linear ranking
%=========================================================================
   
   targ=zeros(n,2);
   
      
           for nn9=1:n
           targ(nn9,:)=[target(nn9,ni,popo) nn9];
           end
           SORTED=sortrows(targ,-1);           
       
       
       POST=zeros(n,P);
   for egyed5=1:n
       for egyed6=1:n
       if SORTED(egyed5,2)==egyed6
       POST(egyed6)=egyed5;
       end
       end
   end
   
   ssp=2;
   for nn4=1:n
        fitness(nn4)=(2-ssp+2*(ssp-1)*((POST(nn4))-1)/(n-1))^2;
   end
%    disp('Fitness es target:')
%    fitness(:,ni,:)
%    target(:,ni,:)
%=========================================================================


     for nn5=1:n  

        if isnan(fitness(nn5))==1
            fitness(nn5)=0.01;       
        end
        if fitness(nn5)<=0
         fitness(nn5)=0.01;
        end
     end

%    
 catch Me
       fclose('all');
       disp(Me)
       disp('Hiba tortent a celfuggveny szamitasakor!')
end


    
 %----------------------------------------------------------------------------------------------------
     MAXIM=max(fitness);
         if 1<MAXIM && MAXIM<=10
            c=10;
         end
         if 0.1<MAXIM && MAXIM<=1
            c=100;
         end
         if MAXIM<=0.1
            c=10000;
         end 
%c
%MAXIM
%p(:,:,ni,2)
   if popo==1
    disp('Szelekcio.')
   end

   % Roulette wheel selection
%------------------------------------------------------------------------------  
   psel=zeros(m,n);
   sel=zeros(n,1);
   % p(:,:,ni,2)
   
%constructing the wheel
       

    fitc=c*(fitness);
     
    if fitc(:)==0
        fitc(:)=0.01;
    end
    rws=sum(round(fitc+1),1);
    RW=zeros(rws,1);
    
   j=round(fitc+1)';
    
       clear l
    for l=1:n   
    if l==1
        
     for  kw=1:j(1);
         RW(kw)=1;
     end 
      kwch=kw;
    else
     for  kw=kwch+1:kwch+j(l);
        RW(kw)=l;
     end
     kwch=kw;
    end
    end
    %RW
    for e=1:n
        sel(e)=randi(rws);
        psel(:,e)=p(:,RW(sel(e))); 
    end


    p=psel;
%---------------------------------------------------------------------------------------------

 %crossover and mutation
%--------------------------------------------------------------------------------------------- 

 % initializing crossover and mutation effect
     

     if mutrate1*m*n<=thresholdM
        mutrate1=thresholdM/(m*n);
     end
     MUT(ni)=mutrate1;
     nmut=round(m*n*mutrate1); 
     sign=[1;-1];
     accur=20;
     Mutshrink=1;
      % probabilistic determination of mutation effect 
     crossover=[ones(round(n*pcr),1);zeros(n-round(n*pcr),1)]; % crossover probability is set to pcr%
     % mutation probability is set to mutrate
     
     Mutmask=zeros(m*n,1);
     for maskelement=1:nmut
         Mutmask(maskelement)=sign(randi(2)); 
     end
   
            
%p(:,:,ni,2)
 %-----------------------------------------------------------------------
   if popo==1
    disp('Rekombinacio.')
   end
 

  if d1<=thresholdD
     d1=thresholdD;
  end
 DT(ni)=d1;
 

% intermediate recombination
%%-------------------------------------------------------------------
tolerance=0;
offspring=zeros(m,ngy);
spsp=1;

while spsp~=n+1
any=randi(n);
ap=randi(n-1);
pany=p(:,any);
%spsp
papsel=cat(2,p(:,1:any-1),p(:,any+1:n));
pap=papsel(:,ap);
%pany,pap
 if crossover(randi(n))==1
 for nno=1:ngy
    for ijj=1:m
            a=(0.5+0.5*randn(1,1)/3)*(-d1)+(0.5+0.5*randn(1,1)/3)*(1+d1);
            %a
            offspring(ijj,nno)=pany(ijj)*a+pap(ijj)*(1-a);
         
            %apply confidence intervals
            if pmin(ijj)>=0
            if offspring(ijj,nno)<pmin(ijj)-tolerance;
            offspring(ijj,nno)=pmin(ijj)-tolerance;
            end
            else
            if offspring(ijj,nno)<pmin(ijj)-tolerance;
            offspring(ijj,nno)=pmin(ijj)-tolerance;
            end            
            end
            if pmax(ijj)>=0
            if offspring(ijj,nno)>tolerance+pmax(ijj);
            offspring(ijj,nno)=tolerance+pmax(ijj);
            end
            else
            if offspring(ijj,nno)>pmax(ijj)+tolerance;
            offspring(ijj,nno)=pmax(ijj)+tolerance;
            end
            end
                   %Zwingen % Constraints: - esetleg k???nyszerfelt???telek (k???s???bb)    
                if flag(ijj)=='i'
                    if (offspring(ijj,nno)-round(offspring(ijj,nno)))<-0.05
                    offspring(ijj,nno)=floor(offspring(ijj,nno));
                    elseif (offspring(ijj,nno)-round(offspring(ijj,nno)))>0.05
                    offspring(ijj,nno)=ceil(offspring(ijj,nno));
                    else
                    offspring(ijj,nno)=round(offspring(ijj,nno));
                    end
                end            
    end

 end
 p(:,spsp:spsp+1)=offspring;
 else 
 p(:,spsp)=pany;
 p(:,spsp+1)=pap;

 end
spsp=spsp+2;
end
if popo==1
pref=p;
end
%%--------------------------------------------------------------------

 %----------------------------------------------------------------------------------
   if popo==1
      disp('Mutacio.')
   end
   
  delt=zeros(accur,1);
  
 %p(:,:,ni+1,2) 
  
  for nnn=1:n
      for mut=1:m
      alfaloc=zeros(accur,1);
      alfaloc(randi(accur))=1;
      for mutvector=1:accur 
      delt(mutvector)=alfaloc(mutvector)*1.4^(-mutvector);
      end
      if pmin(mut)<0
         sign=-1;
      else
         sign=1;
      end
      Delta=sum(delt,1);
      if Delta<0.01
          Delta=0.01;
      end
      %Delta
      pmut=p(mut,nnn);
      eff=Mutmask(randi(m*n))*(pmax(mut)-sign*pmin(mut))*Mutshrink*Delta;
                if flag(mut)=='i'
                    if (eff-round(eff))<-0.05
                    eff=floor(eff);
                    elseif (eff-round(eff))>0.05
                    eff=ceil(eff);
                    else
                    eff=round(eff);
                    end
                end      

        pmut=pmut+eff;
            %apply confidence intervals
            if pmin(mut)>=0
            if pmut<pmin(mut)-tolerance;
            pmut=pmin(mut)-tolerance;
            end
            else
            if pmut<pmin(mut)-tolerance;
            pmut=pmin(mut)-tolerance;
            end            
            end
            if pmax(mut)>=0
            if pmut>tolerance+pmax(mut);
            pmut=tolerance+pmax(mut);
            end
            else
            if pmut>pmax(mut)+tolerance;
            pmut=pmax(mut)+tolerance;
            end
            end
      
       p(mut,nnn)=pmut;
      end
  end
  
% if popo==1
% disp(p-pref)
% end


%% Reinsertion
%%====================================================================
% Random reinsertion
   if popo==1
      disp('Reinzercio')
   end

%p(:,:,ni+1,2)
ptemp=zeros(m,n);

        shuffle=randperm(n);
        for nn8=1:n
        ptemp(:,nn8)=p(:,shuffle(1,nn8));
        end
        p=ptemp;

%p(:,:,ni+1,2)
%--------------------------------------------------------------------------------------------------------------------
p0(:,:,popo)=p;
end

   MAXT=min(min(target,[],1),[],3);
   MEANT=mean(mean(target,1),3);
   
   POPT=sort(reshape(target(:,ni,:),1,n*P));
   %POPNI(:,ni)=POPT(1:sPOP);
   
   PlotPop(ni)=max(POPT)-min(POPT);
   u1=0;  

%----------------------------------------------------------------------------------------
%Break-criteria
%----------------------------------------------------------------------------------------

  %target(:,ni,:)
  %disp('Parameterek:')
    
%clear popo
 
	%p(:,:,ni,2)
%        for egyed2=1:n
%            for popo=1:P
%             if target(egyed2,ni,popo)==MAXT(ni)
%              disp(pdisp(:,:,egyed2,popo))
%              disp(target(egyed2,ni,popo))
%             end
%            end
%        end
	clear egyed2
    %clear popo
    %jjt=1;
    egyed=1;
    popo=1;
    if mod(ni,save_interval)==0 && ni>=2 % save control result in every 10th step
      while egyed<=n
            while popo<=P
            if target(egyed,ni,popo)==MAXT(ni)
             disp(target(egyed,ni,popo))
             %if target(egyed,ni,popo)==POPT(jjt) 
             eredmeny_control(:,:,ni/save_interval,kis)=pdisp(:,:,egyed,popo); 
             %disp(eredmeny_control(:,:,jjt,ni/save_interval,kis))
             save(['eredmeny_control_oc_l_total_',num2str(saveAs),'.mat'],'eredmeny_control');
             save(['celfuggveny_l_total_',num2str(saveAs),'.mat'],'MAXT')  
             %save(['celfuggveny_ls_total_',num2str(saveAs),'.mat'],'POPNI')
             u1=1;
             %jjt=jjt+1;
             egyed=0;
             popo=1;
             break
             %disp(jjt)
            else
             popo=popo+1;
            end               
            end
            
            egyed=egyed+1;
            popo=1;
            
              if u1==1;
              disp('Kontrol eredmeny mentese!')
              break 
              end
      end           
    end
    
%disp('Parameterek:')
%p(:,:,ni,2)
    %size(p)
    
    %Shutdown in case of reaching criterium
    if min(min(target(:,ni,:),[],1),[],3)<=brcr
       for egyed2=1:n
           for popo=1:P
            if target(egyed2,ni,popo)==MAXT(ni)%POPT(jjt)
             eredmeny_control(:,:,ni/save_interval,kis)=pdisp(:,:,egyed,popo); 
             %if jjt==sPOP
             save(['celfuggveny_l_total_',num2str(saveAs),'.mat'],'MAXT')  
             save(['celfuggveny_ls_total_',num2str(saveAs),'.mat'],'POPNI')
             break
             %end
             %jjt=jjt+1;
            end  
            if target(egyed2,ni,popo)==MAXT(ni)
             disp(target(egyed2,ni,popo))
             %u2=1;
             %break
            end
           end
     
%         if u2==1
%             break 
%         end
       end
        disp('Elerte a konvergencia-kriteriumot:')
        disp(ni)
%        matlabpool close
        poolobj=gcp('nocreate');
        delete(poolobj)
        break   
    end

 
%clear popo
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
      
     p0(:,mige(migrant),POPE)=p0(:,migi(migrant),POPI);
     %p(:,migi(migrant),ni+1,POPI)=zeros(m,1,1,1);
     end  
     
 end
end
%p(:,:,ni+1,2)

%migration in ring topology
if strcmp(migtopo,'ring')==1
    migrant1=zeros(migrants,1);
if mod(ni,migtime)==0
    disp('Migracio.')
 
 for popo=1:P-1
     for migrant=1:migrants
      mige(migrant)=randi(n);
      migi(migrant)=randi(n);
     
     p0(:,mige(migrant),popo+1)=p0(:,migi(migrant),popo);
   
     if popo==1
         migrant1(migrant)=mige(migrant);
     end
     end
 end   
 

 for migrant_b=1:migrants
 migrant2=randi(n);
 p0(:,migrant1(migrant_b),1)=p0(:,migrant2,P);
 
 end
end
end
p0(:,:,2)
     if mod(ni,10) == 0
         P_1 = p0;
         save('P_loop_total.mat','P_1')
         %PP0 = PP0 + 1;
     end
 
     %draw fitness
     if mod(ni,5)==0
     plast5(:,:,popo)=p;
     
     hold on
     figure(kis), clf; 
     
     %subplot(3,1,1);
     errorbar(1:ni,MAXT(1:ni)'+PlotPop(1:ni)/2,PlotPop(1:ni)/2,'.')
     hold on,
     scatter(1:ni,MEANT(1:ni),'xr')
     
     ylabel 'Elteres normaja'
     xlabel 'Iteracios lepes'
%      subplot(3,1,2);
%      scatter(1:nii,DT,'ok')
%      ylabel 'd'
%      xlabel 'Iteracios lepes'
%      subplot(3,1,3);
%      scatter(1:nii,MUT,'ok')
%      ylabel 'mutrate'
%      xlabel 'Iteracios lepes'
     title 'Aramhurkok keresese' 
     end
     drawnow;
disp('Lepesszam:')
disp(ni)


end
%==========================output structure==================================== 
disp('Az eredmenyek:')
eredmeny2=eredmeny_control;
disp(eredmeny2)
save('eredmeny.mat','eredmeny2')
%save('eredmeny2.txt','eredmeny2','ascii')
%load('eredmeny.mat');
%matlabpool close

poolobj=gcp('nocreate');
delete(poolobj);

toc
%*********************************************************************************
end


