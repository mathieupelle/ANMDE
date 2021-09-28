function D = get_FourierDiffMatrix(N,NST)
%Create the Fourier differentiation matrix for a given set of nodes
%Use negative sum trick if required)
D = zeros(N,N);

for i = 1:N
    for j = 1:N
        if i == j
            D(i,j) = 0; %The diagonal is zero
        else
            D(j,i) = 0.5*(-1)^(i+j)*cot(pi/(N)*(j-i));
        end
    end
    if NST
        D(j,j) = -sum(D(i,:)); %Diagonal is not zero if NST
    end
end

end