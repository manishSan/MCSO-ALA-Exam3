function [U, B, V] = Bidiag_Francis_Step_Update_U_V(U, B, V)
%Francis_Step Perform one Francis Implicit QR Step

    [ m, n ] = size( B );
    if m < 2 
        return
    end
    
    % we have the form -> Ut * A * V = B
    %                     A = U * B 

    % Introduce the bulge
    % Compute the first Givens' rotation
    G = Givens_rotation( [ B(1,1)^2 - (B(m-1, m)^2 + B( m,m )^2) 
                                B(1,1) * B(1,2) ]);

    % Compute the updates to the matrix
    E_v = eye(m);
    B( 1:2, 1:2 ) = B( 1:2, 1:2 ) * G;
    E_v( 1:2, 1:2 ) = G;
    % (G' * V') = (V * G)'
%     eet_out = E_v*E_v'
    V = V * E_v; 
%     vvt_out = V*V'
%     Vt = V'; 
%     Vt( 1:2, 1:2 ) = G' * Vt( 1:2, 1:2 );
%     V = Vt';

    % Chase the bulge until it is in the last row of the matrix
    for i=1:m-2
        % Applying Givens_rotation from left
        F = Givens_rotation( [ B(i+0, i)
                               B(i+1, i) ]);
 
        B( i:i+2, i:i+2 ) = [F'     [0
                                     0]
                            [0 0]    1] * B( i:i+2, i:i+2 );

        E_u = eye(m);
        E_u(i:i+1, i:i+1) = F;
        eueuT_inside = E_u*E_u'
        U = U * E_u;
        uut_inside = U * U'
%         U( i:i+2, i:i+2 ) = U( i:i+2, i:i+2 ) * ...
%                                     [F      [0
%                                              0]
%                                     [0 0]    1];

        % Applying Givens_rotation from right
        G = Givens_rotation( [ B(i+0, i+1)
                               B(i+0, i+2) ]);
        
        B( i:i+2, i:i+2 ) = B( i:i+2, i:i+2 ) * [ 1    [ 0 0 ]
                                                [ 0
                                                  0 ]     G  ];
        
%         Vt = V'; 
%         Vt( i:i+2, i:i+2 ) = [ 1    [ 0 0 ]
%                              [ 0
%                                0 ]     G' ] * Vt( i:i+2, i:i+2 );
%         V = Vt';
        E_v = eye(m);
        E_v( i+1:i+2, i+1:i+2 ) = G;
%         eet_inside = E_v*E_v'
%         V( i:i+2, i:i+2 ) = V( i:i+2, i:i+2 ) * [ 1    [ 0 0 ]
%                                                 [ 0
%                                                   0 ]     G ];
        V = V * E_v;
%         vvt_inside = V*V'
    end 

    % Remove the bulge from the last row
    L = Givens_rotation( [ B(m-1, m-1)
                           B(m, m-1) ]);

    B( m-1:m, m-1:m ) = L' * B( m-1:m, m-1:m );
    E_u = eye(m);
    E_u(m-1:m, m-1:m) = L;
    eueuT_outside = E_u*E_u'
    U = U * E_u;
    uut_inside = U * U'
%     U( m-1:m, m-1:m ) = U( m-1:m, m-1:m ) * L;
    
end