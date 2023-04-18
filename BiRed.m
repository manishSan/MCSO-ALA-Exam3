function [ A_out, t_out, r_out ] = BiRed( A, t, r )

  [ ATL, ATR, ...
    ABL, ABR ] = FLA_Part_2x2( A, ...
                               0, 0, 'FLA_TL' );

  [ tT, ...
    tB ] = FLA_Part_2x1( t, ...
                         0, 'FLA_TOP' );

  [ rT, ...
    rB ] = FLA_Part_2x1( r, ...
                         0, 'FLA_TOP' );
                     
  while ( size( ATL, 1 ) < size( A, 1 ) )

    [ A00,  a01,     A02,  ...
      a10t, alpha11, a12t, ...
      A20,  a21,     A22 ] = FLA_Repart_2x2_to_3x3( ATL, ATR, ...
                                                    ABL, ABR, ...
                                                    1, 1, 'FLA_BR' );

    [ t0, ...
      tau1, ...
      t2 ] = FLA_Repart_2x1_to_3x1( tT, ...
                                    tB, ...
                                    1, 'FLA_BOTTOM' );

    [ r0, ...
      rho1, ...
      r2 ] = FLA_Repart_2x1_to_3x1( rT, ...
                                    rB, ...
                                    1, 'FLA_BOTTOM' );
                                
    %------------------------------------------------------------%
    % compute householder transform for first column, 
    % and store the householder vector as a21
    % implicitly a21 in the updated matrix equals to zero vector.
    [r1, u21, tau1] = Housev(alpha11, a21);
    alpha11 = r1;
    a21 = u21;
    
    % From section 3.3.4
    % update (a12t / A22) = H((1 / u21), tau1)(a12t / A22)
    % where H = I - 1/tau (u * ut)
    % let w12t = (a12t + u21t * A22) / tau1
    % a12t = a12t - w12t
    % A22 = A22 - u21 * w12t
    w12t = (a12t + (u21') * A22) / tau1;
    a12t = a12t - w12t;
    A22 = A22 - u21 * w12t;

    if (size( a12t', 1 ) > 0)
        % Next introduce zeros in the first row
        % store u12t in a12t
        [u12, rho1] = Housev1(a12t');
        a12t = u12';

        % set the first entry of u12 explicitly to 1
        u12(1) = 1;
    
        % next we update A22 = A22 * H(u12, rho1)
        % A22 := A22 - A22 (u12 * u12')/rho1
    
        v12 = (u12 * (u12'))/rho1;
        A22 = A22 - A22 * v12;
    else 
        rho1 = 0;
    end 
    
    %------------------------------------------------------------%

    [ ATL, ATR, ...
      ABL, ABR ] = FLA_Cont_with_3x3_to_2x2( A00,  a01,     A02,  ...
                                             a10t, alpha11, a12t, ...
                                             A20,  a21,     A22, ...
                                             'FLA_TL' );

    [ tT, ...
      tB ] = FLA_Cont_with_3x1_to_2x1( t0, ...
                                       tau1, ...
                                       t2, ...
                                       'FLA_TOP' );

    [ rT, ...
      rB ] = FLA_Cont_with_3x1_to_2x1( r0, ...
                                       rho1, ...
                                       r2, ...
                                       'FLA_TOP' );
                                   
  end

  A_out = [ ATL, ATR
            ABL, ABR ];

  t_out = [ tT
            tB ];
        
  r_out = [ rT
            rB ];

return