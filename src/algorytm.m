% Wykonuje rozkład Cholesky'ego-Banachewicza zoptymalizowany pamięciowo
% dla macierzy trójdiagonalnej.
% 
% PARAMETRY WEJŚCIOWE:
% A_diag: Wektor zawierający główną przekątną macierzy A.
% A_subdiag: Wektor zawierający elementy podprzekątnej macierzy A (nadprzekątna jest taka sama).
% 
% PARAMETRY WYJŚCIOWE:
% L_diag: Wektor zawierający główną przekątną macierzy L.
% L_subdiag: Wektor zawierający elementy podprzekątnej macierzy L.
% 
% OPIS:
% Funkcja oblicza rozkład Cholesky'ego-Banachewicza macierzy trójdiagonalnej A, zapisując
% wyniki w zoptymalizowanej strukturze pamięciowej. Wynikiem jest macierz L,
% której główna przekątna i podprzekątna są zwracane jako osobne wektory.

function [L_diag, L_subdiag] = algorytm(A_diag, A_subdiag)
    n = length(A_diag);
    L_diag = zeros(n, 1); 
    L_subdiag = zeros(n-1, 1);

    for i = 1:n
        if i == 1
            L_diag(i) = sqrt(A_diag(i));
        else
            L_subdiag(i-1) = A_subdiag(i-1) / L_diag(i-1);
            L_diag(i) = sqrt(A_diag(i) - L_subdiag(i-1)^2);
        end
    end
end