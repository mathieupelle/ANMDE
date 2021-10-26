function ux = LTM_a(N,eps)


%Construct A matrix
A = zeros(N+1,N+1);
for i=1:N+1 %Rows loop
    n = i-1;
    for j=i:N+1 %Columns loop
        p = j; %To easily apply the equations
        u1 = 0; u2 = 0; %Initialize value of the coefficient
        %Calculate coefficient value corresponding to first derivative
        if rem(p+n,2) ~= 0
            u1 = 2*(2*n+1);
        end
        try
        A(i,p+1) = A(i,p+1) + u1;
        end
        %Second derivative
        p = p+1;
        if rem(p+n,2) == 0 
            u2 = 4*(n+0.5)*eps*(p*(p+1)-n*(n+1));
        end
        try
        A(i,p+1) = A(i,p+1)+u2; %Construct element
        end
    end
end

%Get Legendre polynomials
%Evaluating Jacobi polynomials with corresponding alpha, beta and order
X = -1:0.01:1; %Legende polynomials domain
Leg = zeros(N+1,length(X));
for x=1:length(X)
    Leg(:,x) = JacobiP(X(x),0,0,N+1)';
end

%Plug in boundary conditions
for i=1:N+1
    A(end-1,i) = Leg(i,1); % At x=-1
    A(end,i) = Leg(i,end); % At x=1
end


%Construct RHS
f = zeros(N+1,1);
f(1) = -1;

%Solve the system
u = A\f;

%Get the solution 
ux = Leg'*u;