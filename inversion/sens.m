function S=sens(looppar,Nmax,fix,w_d,loop_nr,subname)
%% Function to calculate sensitivity-matrix (Jacobian)
%Basic parameters 
pert=1e-4; %Small perturbation
m = length(looppar); % Nr. of pars to estimate
dp=zeros(m,1);

for j=1:m
    dp(j)=abs(looppar(j)*pert);
end

%Preallocating Jacobian
rows = 2*Nmax*Nmax; %size of [gnm,hnm] SH spectrum
S=zeros(rows,m);

w_p = looppar;
  
%Loop calculating S element-wise
for j=1:m
    pnew = looppar;
    pnew2 = looppar;
    pnew(j)=looppar(j)+dp(j); 
    pnew2(j)=looppar(j)-dp(j);

    % %Sliced parameters for the forward calc. functions
    % %Calling forwards using the unchanged and changed params - a
    [Fmvj, tmp]=subname(pnew,Nmax,loop_nr,fix);
    Fmvj=reshape(Fmvj,[],1);
    [Fmvvj, tmp]=subname(pnew2,Nmax,loop_nr,fix);
    Fmvvj=reshape(Fmvvj,[],1);

    %Calculating S columnwise
    dp(j) = dp(j)/w_p(j);
    delta_d = (Fmvj-Fmvvj)./w_d;
    S(:,j)=delta_d/(2*dp(j));
end

 
end


