function [ B ] = Implicit_bidiag_QR( B )
    [ m, n ] = size( B );
     
    curm = m;   % current last row and column in active matrix
    while curm > 1
        % Check if we need to deflate the last row and column
        if abs( B( curm-1,curm ) ) < ... 
                1.0e-10 * ( abs( B( curm-1,curm-1) ) + ...
                            abs( B( curm, curm ) ) )
             B( curm-1, curm ) = 0;
             curm = curm-1;
             % go to next iteration
             continue  
        end
        
        % determine first zero on superdiagonal, starting from last
        % row/column
        icur = curm-1;
        while icur > 1 
            if abs( B( icur-1,icur ) ) < ... 
                    1.0e-10 * ( abs( B( icur-1,icur-1) ) + ...
                                abs( B( icur, icur ) ) )
                % zero found
                B( icur-1, icur ) = 0;
                break;
            else
                % zero not found
                icur = icur - 1;
            end
        end
        
        [ B(icur:curm,icur:curm) ] = ...
            Bidiag_Francis_Step( B( icur:curm, icur:curm ) );
    end
end