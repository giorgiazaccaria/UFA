function [Rx,tSx,A,dPsi,VQ,V,C] = RandGenSx(J, Q, eSx, flc)
c=2;
while c ~= 1
    [Rx,tSx,A,dPsi,VQ,V,C] = RandGenUltrCor(J, Q, eSx, flc);
    c = isCorrelationMatrix(tSx);  
end

tSx = setNegativeToZero(tSx);

end

% ------------- Local Subfunctions-----------------------------------------

function [Rx,tSx,A,dPsi,VQ,V,C] = RandGenUltrCor(J, Q, eSx, flc)


    % Generate V_Q by randPU
    VQ = randPU(J,Q);

    % If at least 2 variables in the classes of V_Q (flc=1)
    sv = sum(VQ, 1);
       
    if flc==1 
        while any(sv < 2) % At least 2 MVs load onto a factor
            [~, svmax] = max(sv);
            [~, svmin] = min(sv);
            [~, i] = max(VQ(:, svmax));
            VQ(i, svmin) = 1;
            VQ(i, svmax) = 0;
            sv = sum(VQ, 1);
        end
    end

    % generate a random J-tree starting from VQ and including (Q-1) random 
    % agglomerations of subsets are generated 
    V=VQ;
    clag=1:Q;
    for q=1:Q-2
        inj=randi(Q);
        inl=randi(Q);
        while clag(inj)==clag(inl)
            inl=randi(Q);
        end
        V=[V, V(:,clag(inj))+V(:,clag(inl))];
        
        imj=clag==clag(inj);
        iml=clag==clag(inl);
        clag(imj)=Q+q;
        clag(iml)=Q+q;
    end
    V=[V,ones(J,1)];

    % Generate C:
    
    q=0;
    id=0.90;
    a=zeros(2*Q-1,1);
    for q=1:Q
        id=id-1./(2*(Q+2));
        a(q) = id+randn*0.01;
    end
    id=id-0.15;
    for q=1:Q-1
        id=id-1./(2*(Q+2));
        a(Q+q) = id+randn*0.01;
    end

    a(sv==1)=1;
    C=diag(a);

    % Generate A
    A=V*C;

    % Generate Psi
    Rx=A*A';
    dPsi=ones(J,1)-diag(Rx);


    % Generate Rx
    Rx= Rx + diag(dPsi);

    % Generate Sx
    %E = randn(J); %error random normal 
    flge = 1;
    % Generate the Error Matrix: it must be positive semidefinite
    % (correlation matrix) but it can be negative values (otherwise we only
    % add correlation to Rtg increasing the correlation coefficients).
    while flge > 0
        %E = (-ones(J)).^(randi(10,J))*unifrnd(zeros(J),e*ones(J));
        E = unifrnd(zeros(J),eSx*ones(J));
        E = E-tril(E)+triu(E,1)'+eye(J);
        if det(E) >= 0
            flge = 0;
        end
    end

    %E = (E+E')./2; 
    %dE=diag(E);
    %E=E-diag(dE);
    %tSx=Rx+E*eSx;

    tSx=Rx+E - eye(J);
    % Check the positive semidefiniteness of Rn
    eigR = eigs(tSx,J);
    eigR = sort(eigR,'descend');
    if eigR(J) < 0
        tSx = abs(tSx +(0.0001+abs(eigR(J)))*eye(J));
    end
    tSx = abs((diag(diag(tSx))^-(1/2))*tSx*(diag(diag(tSx))^-(1/2)));
    %ldS=log(det(S));

%heatmap(tSx)

end

function isCorrelationMatrix = isCorrelationMatrix(matrix)
    % Check if the matrix is square
    [m, n] = size(matrix);
    if m ~= n
        error('Input matrix must be square.');
    end

    % Check if the matrix is symmetric
    isSymmetric = isequal(matrix, matrix');

    % Check if diagonal elements are all 1
    isDiagonalOnes = all(diag(matrix) == 1);

    % Check if off-diagonal elements are between -1 and 1
    isOffDiagonalInRange = all(all(abs(matrix - eye(n)) <= 1));
    isOffDiagonalInRange1 = all(all(abs(matrix - eye(n)) >=0));

    % Check all conditions
    isCorrelationMatrix = isSymmetric && isDiagonalOnes && isOffDiagonalInRange && isOffDiagonalInRange1 ;
    
end

function resultMatrix = setNegativeToZero(inputMatrix)
    % Ensure input is a matrix
    if ~ismatrix(inputMatrix)
        error('Input must be a matrix.');
    end

    % Set negative elements to zero
    resultMatrix = max(inputMatrix, 0);
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
