function [ rMat ] = vectorRotationFromAtoB( A, B )
    %Written by Gavin Taylor, 2017. MIT License

    %calculate rotation matrix to align b with a - opposite to name

    if size(A,2) == 1
        A = A';
    end
    
    if size(B,2) == 1
        B = B';
    end
    
    %does not require unit inputs
    A = A/sqrt(A(1)^2+A(2)^2+A(3)^2);
    B = B/sqrt(B(1)^2+B(2)^2+B(3)^2);
    
    function ssc = calcSSC(v)
        %calculates skew symetric cross product
        ssc = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
    end
        
    if A == B
        rMat = eye(3);
    else if A == B*-1
        rMat = eye(3)*-1;
        else
            v = cross(A,B);
            rMat = eye(3) + calcSSC(v) + ...
                calcSSC(v)^2*(1-dot(A,B))/(norm(cross(A,B))^2);
    end; end
end



