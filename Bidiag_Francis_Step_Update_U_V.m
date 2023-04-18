function [U, B, V] = Bidiag_Francis_Step_Update_U_V(U, B, V)
%Francis_Step Perform one Francis Implicit QR Step

    [ m, n ] = size( B );
%     if m <= 2 
%         return
%     end
    
    % we have the form -> Ut * A * V = B

    % Introduce the bulge
    % Compute the first Givens' rotation
    G = Givens_rotation( [ B(1,1)^2 - (B(m-1, m)^2 + B( m,m )^2) 
                                B(1,1) * B(1,2) ]);

    % Compute the updates to the matrix

    B( 1:2, 1:2 ) = B( 1:2, 1:2 ) * G;
    V( 1:2, 1:2 ) = V( 1:2, 1:2 ) * G;

    % Chase the bulge until it is in the last row of the matrix
    for i=1:m-2
        % Applying Givens_rotation from left
        F = Givens_rotation( [ B(i+0, i)
                               B(i+1, i) ]);
 
        B( i:i+2, i:i+2 ) = [F'     [0
                                     0]
                            [0 0]    1] * B( i:i+2, i:i+2 );

        Ut = U';
        Ut( i:i+2, i:i+2 ) = [F'     [0
                                     0]
                            [0 0]    1] * Ut( i:i+2, i:i+2 );
        U = Ut';

        % Applying Givens_rotation from right
        G = Givens_rotation( [ B(i+0, i+1)
                               B(i+0, i+2) ]);
        
        B( i:i+2, i:i+2 ) = B( i:i+2, i:i+2 ) * [ 1    [ 0 0 ]
                                                [ 0
                                                  0 ]     G  ];

        V( i:i+2, i:i+2 ) = V( i:i+2, i:i+2 ) * [ 1    [ 0 0 ]
                                                [ 0
                                                  0 ]     G  ];
    end 

    % Remove the bulge from the last row
    L = Givens_rotation( [ B(m-1, m-1)
                           B(m, m-1) ]);

    B( m-1:m, m-1:m ) = L' * B( m-1:m, m-1:m );
    Ut = U';
    Ut( m-1:m, m-1:m ) = L' * Ut( m-1:m, m-1:m );
    U = Ut';


end