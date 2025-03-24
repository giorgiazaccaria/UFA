function [Jtrufa, Aufa, Psiufa, Yufa, fufa, inufa, AGFI] = UFA(X,Q,varargin)
%
% 
% problem minimize discepancy 
%
% D(A,Psi) = log(det(AA'+Psi)) + trace(((AA'+Psi)^-1)*S)
% subject to
% Jtree = J-tree binary with first Q columns a partition in Q columns binary and row stochastic
% A     = loading matrix of UFA
% C = diagonal matrix with non-negative vlaues
% Rx = AA'+Psi
% Rx is ultrametric

%
% INPUT
% Required Parameters
% X (n x J) data matrix or X (J x J) square correlation matrix
% Q = number of classes of the partition of variables. 
%
% Optional parameters:
%
% 'Stats'    ->    Default value: 'off', print the statistics of the fit of the
%                  model.
% 'Stand'    ->    Default value 'on', standardize variables and therefore compute UFA 
%                  on the correlation matrix.
%                  If 'off' does not standardize variables and therefore
%                  compute UFA on the variance-covariance matrix
% 
% 'Rndst'    ->    an integer values indicating the intital random starts.
%                  Default '50' thus, repeat the anaysis 50 times and retain the
%                  best solution.
% 'MaxIter'  ->    an integer value indicationg the maximum number of
%                  iterations of the algorithm
%                  Default '100'.
% 'ConvToll' ->    an arbitrary samll values indicating the convergence
%                  tollerance of the algorithm, Default '1e-9'.
% 'Constr'   ->    Default Constr=[] i.e., no constraints.
%                  Vector Constr (J x 1) indicates for each variable if the
%                  variable is constrained to be in a fixed class. 
%                  Constr(j)= number of class class (c(j) between 1 and K) 
%                  Constr(j)=0 no constraints for variable j.          
% 'NNL'       ->   Defaulf 'off', loadings are free
%                  If 'on' loadings are constrained to be non negative
% 'Vin'       ->   Defaulf Vin=[], initial matrix of the partition of variables will be random
%                  If Vin is given, the solution will start from that
%                  partition of variables.
% 'SOFA'      ->   If 'y' second order factor analysis has to be performed
%                  (for the real data application), default 'n'

% OUTPUT
% Jtrufa (J x 2Q-1) binary membership matrix indicating the classes of
%                   variables and there aggregations
%                   Jtrurufa(j,h)=1 if variable j belongs to class h,  
%                   Jtrufa(j,h)=0 otherwise.
% Aufa   (J x Q)    loading matrix of UFA
% Psiufa (J x 1)    vector of the errors
% Yufa   (n x Q)    Factor scores (if an nxJ matrix is given in input)
% fufa   (1)        discepancy funtion at convergence   
% inufa  (1)        number of time optimal solution was found
% AGFI   (1)        Adjusted Goodness of Fit Index


% n = number of objects
% J = number of variables
%
% OPTIONS
%
% Set optional parameters
%
% Required parameters: X and Q

% initialization

warning('off');
[n,J]=size(X);
opts.disp=0;

% centrering matrix
%Jc=eye(n)-(1./n)*ones(n);

if nargin < 2
   error('Too few inputs');
end

if ~isempty(X)
    if ~isnumeric(X)
        error('Invalid data matrix');
    end  
    if min(size(X)) == 1
    error(message('Ultrametric Factor Analysis:NotEnoughData'));
    end
else
    error('Empty input data matrix');
end

if ~isempty(Q)
    if isnumeric(Q)
        if Q > J 
              error('The number of latent factors larger that the number of variables');
        end    
    elseif Q < 1 
              error('Invalid number of latent factors');
    end
else
    error('Empty input number of latent factors');
end

% Optional parameters   
pnames = {'Stats' 'Stand' 'Rndst' 'MaxIter' 'ConvToll' 'Constr' 'SOFA' 'NNL' 'Vin'};
dflts =  { 'off'    'on'     50     100       1e-9    zeros(J,1) 'n' 'off' [] };
[Stats,Stand,Rndst,MaxIter,ConvToll,Constr,SOFA,NNL,Vin] = internal.stats.parseArgs(pnames, dflts, varargin{:});

if SOFA == 'y'
    HOConstr = [1 1 0 0 1 1 1 1 0 0 0 0; 0 0 1 1 0 0 0 0 1 1 1 1; zeros(J,1)'; zeros(J,1)'; zeros(J,1)']';
end


% Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(Stats)
    if ischar(Stats)
       StatsNames = {'off', 'on'};
       js = strcmpi(Stats,StatsNames);
           if sum(js) == 0
              error(['Invalid value for the ''Statistics'' parameter: '...
                     'choices are ''on'' or ''off''.']);
           end
       Stats = StatsNames{js}; 
    else  
        error(['Invalid value for the ''Statistics'' parameter: '...
               'choices are ''on'' or ''off''.']);
    end
else 
    error(['Invalid value for the ''Statistics'' parameter: '...
           'choices are ''on'' or ''off''.']);
end
% end statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Standardization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(Stand)
    if ischar(Stand)
       StandNames = {'off', 'on'};
       js = strcmpi(Stand,StandNames);
           if sum(js) == 0
              error(['Invalid value for the ''Standardization'' parameter: '...
                     'choices are ''on'' or ''off''.']);
           end
       Stand = StandNames{js}; 
       switch Stand
           
       case 'off'
           if J~=n
            Xs = X-ones(n,1)*mean(X);  
           end
       case 'on'
           if J~=n
            Xs = zscore(X,1);
           else
            dXs=diag(X).^-0.5;        % covariance matrix
            Xs=diag(dXs)*X*diag(dXs);  % correlation matrix
           end
       end
    else  
        error(['Invalid value for the ''standardization'' parameter: '...
               'choices are ''on'' or ''off''.']);
    end
else 
    error(['Invalid value for the ''standardization'' parameter: '...
           'choices are ''on'' or ''off''.']);
end
% end Standardization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rndst %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(Rndst)  
    if isnumeric(Rndst)
       if (Rndst < 0) || (Rndst > 1000) 
       error('Rndst must be a value in the interval [0,1000]');
       end
    else
       error('Invalid Number of Random Starts');
    end
else
    error('Invalid Number of Random Starts')
end
% end Rndst %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MaxIter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(MaxIter)  
    if isnumeric(MaxIter)
       if (MaxIter < 0) || (MaxIter > 1000) 
       error('MaxIter must be a value in the interval [0,1000]');
       end
    else
       error('Invalid Number of Max Iterations');
    end
else
    error('Invalid Number of Max Iterations')
end
% end MaxIter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ConvToll %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(ConvToll)  
    if isnumeric(ConvToll)
       if (ConvToll < 0) || (ConvToll > 0.1) 
       error('ConvToll must be a value in the interval [0,0.1]');
       end
    else
       error('Invalid Convergence Tollerance');
    end
else
    error('Invalid Convergence Tollerance')
end
% end ConvToll %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(Constr)  
    if isnumeric(Constr)
        for j=1:J
            if (Constr(j) < 0) || (Constr(j) > Q) 
                error('Constr must be a value in the interval [0,Q]');
            end
        end  
    else
       error('Invalid Constraint');
    end
else
    error('Invalid Constraint')
end
% end Constr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NonNegative-Loadings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(NNL)
    if ischar(NNL)
       NNLNames = {'off', 'on'};
       jcl = strcmpi(NNL,NNLNames);
           if sum(jcl) == 0
              error(['Invalid value for the ''NonNegative-loadings'' parameter: '...
                     'choices are ''on'' or ''off''.']);
           end
       NNL = NNLNames{jcl}; 
    else  
        error(['Invalid value for the ''NNL'' parameter: '...
               'choices are ''on'' or ''off''.']);
    end
end

NNLFLG=0;
js=strcmpi(NNL,'on');
if sum(js)==1
       NNLFLG=1;
end

% end NNL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute var-covar matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if J~=n
    S=cov(Xs,1);
else
    S=Xs;
end
dS=diag(S);
ldS=log(det(S));
if J~=n
if rank(Xs)<J 
     fprintf('\n \n Attention! Singular Variance-Covariance Matrix \n')
%    Speig=eig(S);
%    Speig=Speig(find(Speig>0.2));
%    ldS=det(diag(Speig));
     ldS=1;
end
    st=(1./n)*sum(sum(Xs.^2));
else
    st=trace(Xs);
end
JJ=[1:J]';
zJ=zeros(J,1);
VC=eye(Q); 

% Start the algorithm      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for loop=1:Rndst
    flg=1;
    while flg > 0
        if isempty(Vin)
            V=randPU(J,Q);      % initial V is radom
        else
            V=Vin;              % initial V is given 
        end
        % if there are Constraints on the variables the initial V has to
        % satisfy
        for j=1:J
            if Constr(j)>0
                V(j,:)=VC(Constr(j),:);
            end    
        end
        flg=sum(find(sum(V)==0));
    end

    A=zeros(J,Q);
    
    for g=1:Q
            ibCg=V(:,g); 
            JCg=[JJ(ibCg==1)];
            %Sg=S(JCg,JCg);
            if sum(ibCg)>1
                %[a,~]=eigs(S(JCg,JCg),1,'lr',opts);
                %[a,y,itr] = ACP(Xs(:,JCg));
                %[a,itr] = AF1(Xs(:,JCg));
                [a,~] = AF1(S(JCg,JCg),NNLFLG);
                A(:,g)=zJ;
                cg=mean(a);
                a=ones(size(a,1),1)*cg;
                A(JCg,g)=a;
            else
                A(JCg,g)=1;
            end
    end
    A0=A;
    %
    % initial Psi
    %
    AA=A*A';
    dAA=diag(AA);
    Psi = diag(dS-dAA);
    Sx=AA+Psi;
    fmin=inf;
    f0=inf;

    % iteration phase
    for it=1:MaxIter
        % update V and A for the first Q levels
        for j=1:J
            posmin=JJ(V(j,:)==1);
            if Constr(j) == 0
                for g=1:Q                
                 V(j,:)=VC(g,:);
                    if sum(V(:,posmin))>0
                        ibCg=V(:,g);           % new class V
                        ibCpm=V(:,posmin);     % old class V
                        JCg=[JJ(ibCg==1)];
                        JCpm=[JJ(ibCpm==1)];
                        if sum(ibCg)>1
                            %[a,~]=eigs(S(JCg,JCg),1,'lr',opts);
                            %[a,itr] = AF1(Xs(:,JCg));
                            [a,~] = AF1(S(JCg,JCg),NNLFLG);

                            if sum(a)<0
                                a=-a;
                            end
                            A(:,g)=zJ;
                            cg=mean(a);
                            a=ones(size(a,1),1)*cg;
                            A(JCg,g)=a;
                        else
                            A(:,g)=zJ;
                            A(JCg,g)=1;
                        end
                        if sum(ibCpm)>1
                            %[aa,~]=eigs(S(JCpm,JCpm),1,'lr',opts);
                            %[aa,itr] = AF1(Xs(:,JCpm));
                            [aa,~] = AF1(S(JCpm,JCpm),NNLFLG);

                            if sum(aa)<0
                                aa=-aa;
                            end
                            A(:,posmin)=zJ;
                            cg=mean(aa);
                            aa=ones(size(aa,1),1)*cg;
                            A(JCpm,posmin)=aa;
                        else
                            A(:,posmin)=zJ;  
                            A(JCpm,posmin)=1;
                        end
                        AA=A*A';
                        dAA=diag(AA);
                        Psi = diag(dS-dAA);
                        Sx=AA+Psi;
                        fQ=log(det(Sx))-ldS + trace(Sx\S)-J;
                        if fQ < fmin
                            fmin=fQ;
                            posmin=g;
                            A0=A;
                        else
                            A=A0;
                        end
                    end
                end
                
            end
            V(j,:)=VC(posmin,:); 
        end
                    
        AA=A*A';
        dAA=diag(AA);
        Psi = diag(dS-dAA);
        Sx=AA+Psi;
        fQ=log(det(Sx))-ldS + trace(Sx\S)-J;

        sA=sum(A);
        sV=sum(V);
        ck=sA./sV;
        [~,icc]=sort(ck, 'descend');
        A=A(:,icc);
        V=V(:,icc);
        %GFA = A; %for group factor analysis
        
        VQ = V;
        AQ=A;
        [J,Q] = size(V);
        Ir = zeros(Q-1,1);
        Ic = zeros(Q-1,1);
        lfm = zeros(Q-1,1);
        PP = zeros(Q,1);
        
        % Update A and V for the remaining Q-1 levels
        %
        Rb=pinv(VQ)*S*pinv(VQ').*(ones(Q)-eye(Q));
        Rw=pinv((VQ'*VQ)^2-VQ'*VQ)*diag(diag(VQ'*(S-eye(J))*VQ));
        sVQ=sum(VQ);
        if find(sVQ==1) ~= 0
           Rw(sVQ==1,sVQ==1) = 1;
        end
        Rw=diag(diag(Rw));
        
        %Ac= V*Rw;
        Ac=A.^2;
        
        % Class Leaders vector PP
        PP=1:Q;
        Rbt=triu(Rb,1);
        % Update V
        for q = 1:Q-1
            % Max Correlation between groups "lfm"
            lfm(q) = max(max(Rbt));
                
            [ir,ic]=find(round(Rbt,10)==round(lfm(q),10));
            Rbt(ir,ic)=0;
            ir=ir(1);
            ic=ic(1);
        
            Jvh=[JJ(V(:,PP(ir))==1)]; %indices of v_h
            Jvq=[JJ(V(:,PP(ic))==1)]; %indices of c_q

            while PP(ir) == PP(ic)
                Rbt(ir,ic)=0;
                lfm(q) = max(max(Rbt));
                [ir,ic]=find(round(Rbt,10)==round(lfm(q),10));
                ir=ir(1);
                ic=ic(1);
                Jvh=[JJ(V(:,PP(ir))==1)]; %indices of v_h
                Jvq=[JJ(V(:,PP(ic))==1)]; %indices of c_q
            end

            % if there are Constraints on the variables the initial V has to
            % satisfy
            for j=1:J
               if SOFA == 'y'
                   V(j,Q+q)= HOConstr(j,q);
               end    
            end
            vhq=V(:,PP(ir))+V(:,PP(ic));
                        
            if SOFA == 'y'
                for j=1:J
                   V(j,Q+q)=HOConstr(j,q);
                end
            else
                V=[V vhq];
            end

            
            % Update the class leaders
            Ir(q) = PP(ir);
            Ic(q) = PP(ic);
        
            icPP=PP==PP(ic);
            irPP=PP==PP(ir);
            PP(icPP)=Q+q;
            PP(irPP)=Q+q;
        end 
        %
        % Update Ac 
        for q=1:Q-1      
            irQpq=[];
            icQpq=[];
            irQpq=find(Ir==Q+q);
            if isempty(irQpq)
                itQpq=find(Ic==Q+q);
            else
                itQpq=irQpq;
            end
            if Ir(q)<=Q
                Jvh=[JJ(V(:,Ir(q))==1)]; %indices of v_h
                Ac(Jvh,Ir(q))=real((Ac(Jvh,Ir(q))-lfm(q)).^0.5);  
            end
            if Ic(q)<=Q
                Jvq=[JJ(V(:,Ic(q))==1)]; %indices of v_q
                Ac(Jvq,Ic(q))=real((Ac(Jvq,Ic(q))-lfm(q)).^0.5);    
            end
            if  itQpq >= 1
                Ac(:,Q+q)=real(V(:,Q+q)*(lfm(q)-lfm(itQpq)).^0.5); 
            else 
                Ac(:,Q+q)=real(vhq*(lfm(q)).^0.5);      
            end
        end
        %
        %
        sVQ=sum(VQ);
        icVQ=JJ(sVQ==1);
        Ac(:,icVQ)=VQ(:,icVQ);

        A=Ac;
        if SOFA == 'y'
            V(:,2*Q-1) = zeros(J,1);
            A(:,2*Q-1) = zeros(J,1);
        end
        %
        % Check convergence
        AA=A*A';
        dAA=diag(AA);
        Psi = diag(dS-dAA);
        Sx=AA+Psi;
        f=log(det(Sx))-ldS + trace(Sx\S)-J;
        fdif = f0-f;
            if fdif > ConvToll 
                f0=f; fmin=fQ;A=AQ;A0=A;V=VQ; 
            else
                fprintf('UFA: Loop=%g, Discrepancy= %g, number iteration = %g, fdif=%g\n',loop,f, it,fdif)
                break
            end
    end
    
    % Store the optimal solution
    if loop==1
            Jtrufa=V;
            Aufa=A;
            Psiufa=Psi;
            RxU=Sx;
            if J~=n
                Yufa=zscore(Xs*Aufa,1)*(Aufa'*Aufa)^0.5;  
            else
                Yufa=[];
            end
            fufa=f;
            loopdfa=1;
            inufa=it;
            fdifo=fdif;
            RxU=Sx;
   end
   if f < fufa
       Jtrufa=V;
       fufa=f;
       Aufa=A;
       Psiufa=Psi;
       RxU=Sx;
       if J~=n
        Yufa=zscore(Xs*Aufa,1)*(Aufa'*Aufa)^0.5; 
       end
       loopdfa=loop;
       inufa=it;
       fdifo=fdif;
   elseif abs(f-fufa) < ConvToll && loop>1
       break
   end
end

% sort the final solution in descend order of variance
if J~=n
 varYdfa=var(Yufa,1);
 [~,ic]=sort(varYdfa, 'descend');
 Aufa=Aufa(:,ic);
 Jtrufa=Jtrufa(:,ic);
 Yufa=Yufa(:,ic); 
end

% compute the variace of the second component  

e2k=zeros(Q,1);
cro=zeros(Q,1);
for g=1:(2*Q-1)
    ibCg=Jtrufa(:,g); 
    if sum(ibCg)>1
        JCg=[JJ(ibCg==1)];
        cro(g)=CronbachAlpha(Xs(:,JCg));
        Sg=S(JCg,JCg);
        [U,L]=eig(Sg);
        [l,il]=sort(diag(L), 'descend');
        L=diag(l);
        U=U(:,il);
        e2k(g)=L(2,2);
    else
        e2k(g)=0;
        cro(g)=1;
    end
end

% Performance metrics
Sxufa=Aufa*Aufa'+ Psiufa;
invSxS=Sxufa\S;

% For the real data application
%n = 277;
%lldfa=-(1./2)*n*log(2*pi)-(1./2)*n*(log(det(Sxufa))+ trace(invSxS));
%BIC=-2*lldfa+log(n)*pm;
%AIC=-2*lldfa+2*pm;

GFI= 1- trace(((Sxufa^-1)*S-eye(J))^2)/trace((Sxufa^-1*S)^2);
df = J*(J-1)/2 -2*Q +1;
AGFI = 1 - (J*(J-1)/(2*df))*(1-GFI);

% Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

js=strcmpi(Stats,'off');
if sum(js) == 0
    fprintf('\n \n Ultrametric FA (Final) \n Unknown Parameters =%g,  Discrepancy =%g, \n Goodness of Fit =%g \n Adjusted GOF =%g, \n loopdpca=%g, iter=%g, fdif=%g \n \n Chi-square =%g,  p-value=%g \n degree of freedom =%g, \n BIC=%g, \n AIC=%g, \n Root-Mean-Square-Error-of-Approximation \n RMSEA=%g \n',pm,fufa,GFI,AGFI, loopdfa, inufa,fdifo,X2,pValue,df, BIC,AIC, RMSEA)
    fprintf('\n \n Ultrametric Factor Analysis results\n')
    fprintf('\n \nFactor            Variance Explained     Perc. Expl. Var.   Cumulated Var.  Perc. Expl. Var.\n')
    VarCum=0;
    for k=1:Q
        VarCum = VarCum + sum(Aufa(:,k).^2);
        fprintf('factor   (%2g)      %f               %f          %f        %f  \n', k, sum(Aufa(:,k).^2), sum(Aufa(:,k).^2)./st*100, VarCum, VarCum./st*100)
    end
    fprintf('\n \n Loading Matrix (Manifest Variables x Latent Variables)\n')

    for k=1:Q
        fprintf('Unidimensionality Assessment: variance of the second component of class(%g)    %f\n', k,e2k(k))
    end
    for k=1:(2*Q-1)
        fprintf('Reliability Assessment: Alpha of Cronbach of class(%g)    %f\n', k,cro(k))
    end
    com=sum(Aufa'.^2);
    fprintf('\n \n         Path           Corr coef       Std Error        Pr(p>|Z|)        Var Error        Communality\n' )
    for j=1:J
        p=zeros(J,1);
        %p=zeros(J,2);
        for k=1:Q
            %p(j,:)= normcdf([-1*abs(AA(j,1)/SD(j,j).^0.5) abs(AA(j,1)/SD(j,j).^0.5)]);
            if Aufa(j,k)~=0
                p(j)= r2pv(Aufa(j,k),n);
                fprintf('X(%g) <-- factor(%g)     % f        %f         %f         %f         %f\n', j, k, Aufa(j,k), Psiufa(j).^0.5./sqrt(n), p(j), Psiufa(j), com(j))    
            end
        end
    end
end

fprintf('UFA Final: Loop=%g, Discrepancy= %g, Goodness of Fit = %g, Adjusted GOF = %g, iter=%g, fdif=%g\n',loop,f, GFI, AGFI, it,fdif)
end

% ------------- Local Subfunctions-----------------------------------------
function [a,y,itr] = ACP(Xr,NNLFLG)
maxit=300;
[n,Q]=size(Xr);
a = rand(Q,1); 
tol = 1e-9;
error = inf;
last = inf;
itr = 0;
while ~(abs(last-error)<error*tol) && itr<=maxit
    itr = itr+1;

    y = Xr*a./(a'*a);
    a=Xr'*y./(y'*y);
    if NNLFLG == 1   % if loadings must be non negative
        psa=find(a>=0);                                 %
        if size(psa) < Q                                %
            y=Xr(:,psa)*a(psa)./(a(psa)'*a(psa));       %
            as=a<0;                               % if non negative
            a(psa)=Xr(:,psa)'*y./(y'*y);                %
            a(as)=0;                                    %
        end                                             %
    end
        last = error;    
        e = y-Xr*a;
        error = e'*e/n;
    end
a=a./sqrt(a'*a);
y=Xr*a;
end

function [A,itr] = AF1(S, NNLFLG)
%function [A,itr] = AF1(Xr);
% Factor Analysis
%    Sx=AA' + Psi
maxit=300;
dS=diag(S);
Psi=diag(dS); 

% given initial Psi compute A
Psii=diag((dS.^-0.5));
%[U,L]=eigs(Psii*S*Psii,1);
D=Psii*S*Psii;
[U,~,~] = ACP(D, NNLFLG);
L=U'*D*U;
Ps=diag((dS.^0.5));
ll=diag(L);
ll=(ll-1).^0.5;
A=Ps*U*diag(ll);

% compute discrepancy
AA=A*A';
Sx=AA+Psi;
discrepO =log(det(Sx))+trace((Sx^-1)*S);


tol = 1e-10;
discrep = inf;
itr = 0;

while ~(abs(discrep-discrepO)< tol) && itr<= maxit
    itr = itr+1;
    if itr ~=1
        discrepO=discrep;
    end
    % given A update Psi
    AA=A*A';
    dAA=diag(AA);
    Psi = diag(dS-dAA);

    % given Psi update A
    dPsi=diag(Psi);
    Psii=diag((dPsi.^-0.5));
    D=Psii*S*Psii;
    [U,~,~] = ACP(D,NNLFLG);
    %[U,ll]=eigs(D,1);
    L=U'*D*U;
    Ps=diag((dPsi.^0.5));
    ll=diag(L);
    ll=ll-1;
    ll=ll.^0.5;
    A=Ps*U*diag(ll);
    AA=A*A';

    % compute discrepancy
    Sx=AA+Psi;
    discrep=log(det(Sx))+trace((Sx^-1)*S);
end
end


function p=r2pv(r,n)
%
% 	p=r2pv(r,n)
%
% r = estimated correlation coefficient (IE |r| <= 1)
%   = (1/n)*(x'*y) for col vectors x,y of length n
% n = no. samples used
% p = P-value based on |r| (two sided) with rho=0 (null case)
%
% NOTES: following Cramer, p.400, convert r to a t and use what we have for t 
if n < 3
    error('n < 3');
end
if r==1. 
    p=0; 
    return;
end
t=sqrt(n-2)*r/(sqrt(1-r*r)); 	% this is t with n-2 d.f.
t=abs(t);							% use |t| for two sided P-value
p=2*(1-tcdf(t,n-2));
end

function mri=mrand(N)
%
% modified rand index (Hubert & Arabie 1985, JCGS p.198)
%
n=sum(sum(N));
sumi=.5*(sum(sum(N').^2)-n);
sumj=.5*(sum(sum(N).^2)-n);
pb=sumi*sumj/(n*(n-1)/2);
mri=(.5*(sum(sum(N.^2))-n)-pb)/((sumi+sumj)/2-pb);
end

function [as,varargout] = CronbachAlpha(x)
% CronbachAlpha
% 
% Description:	calculate Cronbach's alpha for a set of psychometric measurements
% 
% Syntax:	[as,au] = CronbachAlpha(x)
% 
% In:
% 	x	- an nRep x nItem array of ratings, so that each row is the set of
%		  obvservations from one repetition and each column is the set of all
%		  observations for a given item
% 
% Out:
% 	as	- the standardized Cronbach's alpha
%	au	- the unstandardized Cronbach's alpha
% 
% Updated: 2012-09-24
% Copyright 2012 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
nItem	= size(x,2);

%logical array for selecting upper triangular part of the correlation and
%covariance matrices, where the good stuff is
	b	= triu(true(nItem),1);

%standardized alpha
	%pairwise correlations between items
	if size(x,1) == nItem
    	r = x;
    else
        r	= corrcoef(x);
    end
	%mean of the meaningful, non-redundant correlations
		r	= mean(r(b));
	
	as	= nItem*r/(1 + (nItem-1)*r);

%unstandardized alpha
if nargout>1
	%variance/covariance matrix
		vc	= cov(x);
	%mean variance (variances are along the diagonal)
		v	= mean(diag(vc));
	%mean covariance, not including variances
		c	= mean(vc(b));
	
	varargout{1}	= nItem*c/(v + (nItem-1)*c);
end
end

function [U]=randPU(n,c)

% generates a random partition of n objects in c classes
%
% n = number of objects
% c = number of classes
%
U=zeros(n,c);
U(1:c,:)=eye(c);

U(c+1:n,1)=1;
for i=c+1:n
    U(i,1:c)=U(i,randperm(c));
end
U(:,:)=U(randperm(n),:);
end


