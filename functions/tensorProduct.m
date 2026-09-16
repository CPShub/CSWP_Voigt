function T = tensorProduct(A,B)
    % Check that matrices are square
    [rowA,colA] = size(A);
    [rowB,colB] = size(B);

    if ~((rowA == colA) && (rowB == colB))
        error("Matrices are not square!")
    end

    if rowA ~= rowB
        error("Matrices are not of same shape")
    end

    dim = rowA;

    T = 0.5 * ( ...
          reshape(A,dim,1,dim,1) .* reshape(B,1,dim,1,dim) ...
        + reshape(A,dim,1,1,dim) .* reshape(B,1,dim,dim,1) );
end