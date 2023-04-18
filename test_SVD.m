format long

% Seed the random number generator so we always get the same random numbers
rng(1)

% Set the matrix size
m = 2;

% Generate a random matrix
A = rand( m, m );

% Compute the singular values using matlabs intrinsic function
sigmas = svd( A );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reduce the matrix to bidiagonal form.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create vector in which to store the scalars tau from the Householder
% transformations
t = rand( m, 1 );

% Create vector in which to store the scalars rho from the Householder
% transformations.  (Yes, actually, we only need m-1 scalars...)
r = rand( m, 1 );

% Compute the reduction to tridiagonal form.  You need to write BiRed.  See
% "Ponder This 11.2.3.3
[ B, t, r ] = BiRed( A, t, r );

% Quick check if it was probably done correctly: Check the singular values of
% the bidiagonal matrix (extracted from B) with those of the original matrix.

% Extract the bidiagonal matrix from B.
Bi = BiFromB( B );

disp('Norm of difference between singular values of A and singular values of Bi')
disp( norm( svd( A ) - svd( Bi ) ) );
svd(A)
svd(Bi)
if norm( svd( A ) - svd( Bi ) ) > 1.0e-14 
    disp( 'Something is wrong with your reduction to bidiagonal form' ); 
    disp( 'left colum:   singular values of A.' );
    disp( 'right column: singular values of Bi.');
    [ svd( A )  svd( Bi ) ]
end

% Do this once the other parts are working.
% Form the matrices U_A and V_A that are stored in B below the diagonal and
% above the first superdiagonal so that A = U_A * Bi * V_A.
% Hints: 
% - Use the routine FormQ from Week 3!
% - To form U_A: 
%   Where are the Householder vectors stored?
% - To form V_A:
%   Where and how are the Householder vectors stored?  
%   What do you need to do to put them in the right format to pass to FormQ?  
%   Where in V_A does the resulting unitary matrix go?
% Uncomment the lines

% For U_A, the householder vectors are stored in the column
U_A = FormQ(B, t);

% For V_A, the householder vectors are stored in the rows
% so we take a transpose of B
V_A = eye( m, m );
Bt = B';
V_A( 2:m, 2:m ) = FormQ(Bt( 2:m, 1:m-1), r(1: m-1));

% Uncomment the below to check if you got it right:
disp( 'Difference between U_A * Bi * V_A and A:' );
A
U_A * Bi * V_A'
norm( U_A * Bi * V_A' - A )
if ( norm( U_A * Bi * V_A' - A ) > 1.0e-10 )
    disp( 'Something is wrong with your U_A and V_A' ); 
    disp( 'Difference between U_A * Bi * V_A and A:' );
    U_A * Bi * V_A' - A;  
end

% Write a single step of a "bidiagonal Francis Step" that introduces the
% bulge and chases it out, working with the bidiagonal matrix in Bi.
% I suggest you do this in stages.
% 
% First write it so that you compute the Givens' rotations and you appy
% these to the entire rows and columns to which they are to be applied.  In
% other words, literally do what is in 11.2.4, ignoring where the zeroes
% are when you apply the Givens' rotation.
%
% Next, think carefully about where the nonzeroes are (or are introduced)
% and only update those elements in the rows and/or columns.
%
% Finally, and this is optional, write the routine so you don't corrupt any
% entries below the diagonal and above the first superdiagonal.  Don't
% waste too much time on this.
% Below, we just execute a bunch of these steps, to see if the elements on
% the superdiagonal converge the way they should.

Bi_next = Bidiag_Francis_Step( Bi )
Bi_next = Bidiag_Francis_Step( Bi_next )
Bi_next = Bidiag_Francis_Step( Bi_next )
Bi_next = Bidiag_Francis_Step( Bi_next )
Bi_next = Bidiag_Francis_Step( Bi_next )
Bi_next = Bidiag_Francis_Step( Bi_next )
Bi_next = Bidiag_Francis_Step( Bi_next )
Bi_next = Bidiag_Francis_Step( Bi_next )
Bi_next = Bidiag_Francis_Step( Bi_next )
Bi_next = Bidiag_Francis_Step( Bi_next )
Bi_next = Bidiag_Francis_Step( Bi_next )
Bi_next = Bidiag_Francis_Step( Bi_next )

% We wrote the following routine, which wraps deflation in a loop around
% the routine Bidiag_Francis_Step that you wrote.

D = Implicit_bidiag_QR( Bi );

% Note: in the below, we make sure D is diagonal by extracting its
% diagonal.  We do svd( diag ( diag( D ) ) to order the singular values (on
% the diagonal of D).

disp('Norm of difference between diagonal elements of D and singular values of A')
disp( norm( svd( A ) - svd( diag( diag( D ) ) ) ) );

if norm( svd( A ) - svd( diag( diag( D ) ) ) ) > 1.0e-10 
    disp( 'Something is wrong with your Implicit_bidiag_QR' ); 
    [ svd( A )  svd( diag( diag( D ) ) )  ]
end

% Now, modify Implicit_bidiag_QR to create Implicit_bidiag_QR_SVD that also
% updates U_A and V_A by applying the Givens' rotations appropriately,
% which means in the end you get the SVD of A, but with the singular values
% not ordered in the correct order.

% Uncomment to test

% Call your routine that updates U_A and V_A to get the SVD (modulo the

% ordering of the singular values)
[ U, Sigma, V ] = Implicit_bidiag_QR_SVD( U_A, Bi, V_A );

if norm( svd( A ) - svd( diag( diag( Sigma ) ) ) ) > 1.0e-10 
    disp( 'Something is wrong with your Implicit_bidiag_QR_SVD' ); 
    [ svd( A )  svd( diag( diag( Sigma ) ) )  ]
end

disp( 'Difference between A and U * Sigma * transpose(V) ')
disp( norm( A - U * Sigma * V', 'fro' ) );

A - U * Sigma * V'

U
Sigma
V

[U_f,S_f,V_f] = svd(A)