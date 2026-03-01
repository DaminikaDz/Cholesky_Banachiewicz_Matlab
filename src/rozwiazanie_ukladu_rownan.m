% Rozwiązuje układ równań Ax = b
% przy użyciu rozkładu Cholesky'ego-Banachewicza dla macierzy trójdiagonalnych.
% 
% PARAMETRY WEJŚCIOWE:
% A_diag: Wektor zawierający główną przekątną macierzy A.
% A_subdiag: Wektor zawierający elementy podprzekątnej macierzy A (nadprzekątna jest taka sama).
% b: Wektor prawych stron równania Ax = b.
% 
% PARAMETRY WYJŚCIOWE:
% x: Wektor rozwiązania układu równań Ax = b.
% 
% OPIS:
% Funkcja wykorzystuje zoptymalizowaną pamięciowo wersję rozkładu
% Cholesky'ego-Banachewicza dla macierzy trójdiagonalnych. Obliczenia są prowadzone
% w dwóch krokach: najpierw rozwiązanie L^T y = b,
% następnie rozwiązanie L x = y.
function x = rozwiazanie_ukladu_rownan(A_diag, A_subdiag, b)
    [L_diag, L_subdiag] = algorytm(A_diag, A_subdiag);
    n = length(b);

    % Rozwiązanie L^T y = b
    y = zeros(n, 1);
    y(1) = b(1) / L_diag(1);
    for i = 2:n
        y(i) = (b(i) - L_subdiag(i-1) * y(i-1)) / L_diag(i);
    end

    % Rozwiązanie L x = y
    x = zeros(n, 1);
    x(n) = y(n) / L_diag(n);
    for i = n-1:-1:1
        x(i) = (y(i) - L_subdiag(i) * x(i+1)) / L_diag(i);
    end
end