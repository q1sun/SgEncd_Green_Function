function result = Local_Basis_1D(x,vertices,basis_type,basis_index,derivative_degree)

format short e

% evaluation of the (basis_index)-th local basis functions on coordinate x, where derivative_degree = the derivative degree of the basis function
if basis_type == 101 % linear element
    
    if derivative_degree == 0
        
        if basis_index == 1
            result = (vertices(2)-x) / (vertices(2)-vertices(1));
        elseif basis_index == 2
            result = (x-vertices(1)) / (vertices(2)-vertices(1));
        end

    elseif derivative_degree == 1
        
        if basis_index == 1
            result = 1 / (vertices(1)-vertices(2));
        elseif basis_index == 2
            result = 1 / (vertices(2)-vertices(1));
        end
        
    end
    
elseif basis_type == 102 % quadratic element
    
    denominator = (vertices(1)-vertices(2))^2;
    
    if derivative_degree == 0
        
        if basis_index == 1
            result = 2*x.^2 - (vertices(1)+3*vertices(2))*x + vertices(2)^2 + vertices(1)*vertices(2);
        elseif basis_index == 2
            result = 2*x.^2 - (3*vertices(1)+vertices(2))*x + vertices(1)^2 + vertices(1)*vertices(2);
        elseif basis_index == 3
            result = -4*x.^2 + 4*(vertices(1)+vertices(2))*x - 4*vertices(1)*vertices(2);
        end

    elseif derivative_degree == 1
        
        if basis_index == 1
            result = 4*x - (vertices(1)+3*vertices(2));
        elseif basis_index == 2
            result = 4*x - (3*vertices(1)+vertices(2));
        elseif basis_index == 3
            result = -8*x + 4*(vertices(1)+vertices(2));
        end

    elseif derivative_degree == 2
        
        if basis_index == 1
            result = 4;
        elseif basis_index == 2
            result = 4;
        elseif basis_index == 3
            result = -8;
        end
        
    end
    
    result = result/denominator;
    
end

end


