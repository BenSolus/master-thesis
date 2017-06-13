#!/usr/bin/octave

close all;

A = zeros(4, 4);
B = zeros(4, 4);

A(1, 1) = -8.38212e-02; A(1, 2) = 9.50307e-01; A(1, 3) = 2.33956e-01; A(1, 4) = 1.0;
A(2, 1) = 8.90248e-01;  A(2, 2) = 5.09489e-02; A(2, 3) = 1.0;
A(3, 1) = -7.50751e-02; A(3, 2) = 1.0;
A(4, 1) = 1.0;

B(1, 1) = 1.0;
B(2, 1) = 7.50751e-02;  B(2, 2) = 1.0;
B(3, 1) = -8.94073e-01; B(3, 2) = -5.09489e-02; B(3, 3) = 1.0;
B(4, 1) = 2.21651e-01;  B(4, 2) = -9.38387e-01; B(4, 3) = -2.33956e-01; B(4, 4) = 1.0;

A

B

I = A * B;

I
