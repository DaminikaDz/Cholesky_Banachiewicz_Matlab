% Ten skrypt testuje działanie funkcji rozwiązywania układów równań Ax = b 
% dla macierzy trójdiagonalnych przy użyciu zoptymalizowanego pamięciowo 
% rozkładu Cholesky'ego-Banachewicza(sa przechowywane tylko dwie przekątne 
% macierzy A i L, gdzie A=L^T*L)
% 
% Brak parametrów wejściowych i wyjściowych - wszystkie wyniki są 
% wyświetlane w MATLAB Console w formie tabeli.
%
% OPIS:
% Skrypt testuje 6 ciekawych przypadków dla macierzy 
% trójdiagonalnych i odpowiadających wektorów b. Wyniki rozwiązania 
% numerycznego są porównywane z wbudowaną funkcją MATLABa (A \ b).
% 
% Obliczane są także takie metryki jak: 
% - wyznacznik macierzy A (zarówno metodą rozkładu Cholesky'ego-Banachewicza, jak i MATLAB-em),
% - błąd numeryczny wyznacznika,
% - współczynniki stabilności i poprawności,
% - błąd względny rozwiązania,
% - wskażnik uwarunkowania macierzy,
% - błąd rozkłądu.

matrix = {
    {[4; 3; 5; 4], [1; 0; 2]}, ... % Zwykłe wartości
    {[1; 10; 100; 1000; 1e4; 1e5], [1; 5; 10; 50; 100]}, ... % Rosnące wartości
    {[1e6; 2e6; 3e6; 4e6], [1e-6; 1e-5; 1e-4]}, ... % Duże wartości na przekątnej, pozostałe małe
    {[1e9; 1e8; 1e7; 1e6], [1e5; 1e4; 1e3]}, ... % Duże wartości
    {[1; 1; 1; 1], [1; 1; 1]}, ...                % Nie jest dodatnio okreslona
    {[0.1; 0.2; 0.4; 0.8], [0.05; 0.02; 0.06]}, ... % Małe wartości
    {[1e-5; 1e-4; 1e-3; 1e-6], [1e-4; 1e-5; 1e-6]} ... % Skrajnie małe wartości
    
};

% odpowiednie wektory b 
bs = {[1, 1, 1, 1];
    [0.1, 0.01, 0.05, 0.006, 0.0001, 0.00006];
    [1e-6, 2e-6, 3e-6, 4e-6];
    [1e9, 1e8, 1e7, 1e6];
    [1, 1, 1, 1];
    [0.1, 0.02, 0.01, 0.01];
    [1, 10, 0.01, 0.0001]
};
% Inicjalizacja tabeli wyników
results_table = cell(length(matrix), 9);

for i = 1:length(matrix)
    A_diag = matrix{i}{1};
    A_subdiag = matrix{i}{2};
    b = bs{i};

    % Rozwiązanie przy użyciu własnej implementacji
    x_custom = rozwiazanie_ukladu_rownan(A_diag, A_subdiag, b);

    % Odtworzenie pełnej macierzy A do porównania
    n = length(A_diag);
    A = diag(A_diag) + diag(A_subdiag, 1) + diag(A_subdiag, -1);
    
    % Rozwiązanie przy użyciu funkcji wbudowanej MATLAB
    x_builtin = A \ b';

    % Obliczenie wyznacznika przy użyciu rozkładu Cholesky'ego-Banachewicza
    [L_diag, L_subdiag] = algorytm(A_diag, A_subdiag);
    det_A_cholesky = prod(L_diag)^2; % Wyznacznik z rozkładu Cholesky'ego

    % Oblicz wyznacznik przy użyciu funkcji wbudowanej MATLAB
    det_A_builtin = det(A);

    % Obliczenie błędu wyznacznika
    det_error = abs(det_A_cholesky - det_A_builtin);

    cond_A = cond(A);

    % Błąd dekompozycji
    A_reconstructed = diag(L_diag) * diag(L_diag)' + diag(L_subdiag, -1) * diag(L_subdiag, -1)';
    decomposition_error = norm(A - A_reconstructed) / norm(A);

    % Względny błąd rozwiązania
    relative_error = norm(x_custom - x_builtin) / norm(x_builtin);

    x_exact = x_builtin; 
    stability_coefficient = norm(x_custom - x_exact) / (norm(x_exact) * cond_A);
    correctness_coefficient = norm(b - A * x_custom) / (norm(A) * norm(x_custom));

    results_table{i, 1} = sprintf('[%s]', strjoin(arrayfun(@(x) num2str(x, '%.4f'), x_custom, 'UniformOutput', false), ', '));

    results_table{i, 2} = sprintf('[%s]', strjoin(arrayfun(@(x) num2str(x, '%.4f'), x_builtin, 'UniformOutput', false), ', '));
    
    results_table{i, 3} = det_A_cholesky; 
    results_table{i, 4} = det_A_builtin;
    results_table{i, 5} = det_error; 
    results_table{i, 6} = cond_A; 
    results_table{i, 7} = decomposition_error;
    results_table{i, 8} = relative_error; 
    results_table{i, 9} = stability_coefficient; 
    results_table{i, 10} = correctness_coefficient; 

    fprintf('Przykład %d:\n', i);
    fprintf('Rozwiązanie_Własne: %s\n', mat2str(x_custom, 4));
    fprintf('Rozwiązanie_Wbudowana_f: %s\n', mat2str(x_builtin, 4));
    fprintf('Wyznacznik_Cholesky: %.4e\n', det_A_cholesky);
    fprintf('Wyznacznik_Wbudowana_f: %.4e\n', det_A_builtin);
    fprintf('Błąd_Wyznacznika: %.4e\n', det_error);
    fprintf('Wskaźnik_Uwarunkowania: %.4e\n', cond_A);
    fprintf('Błąd_Dekompozycji: %.4e\n', decomposition_error);
    fprintf('Błąd_Względny_Rozwiązania: %.4e\n', relative_error);
    fprintf('Współczynnik_Stabilności: %.4e\n', stability_coefficient);
    fprintf('Współczynnik_Poprawności: %.4e\n\n', correctness_coefficient);
end

disp('Tabela wyników:');
disp(cell2table(results_table, ...
    'VariableNames', { ...
    'Rozwiązanie_Własne', ...
    'Rozwiązanie_Wbudowane', ...
    'Wyznacznik_Cholesky', ...
    'Wyznacznik_Wbudowany', ...
    'Błąd_Wyznacznika', ...
    'Wskaźnik_Uwarunkowania', ...
    'Błąd_Dekompozycji', ...
    'Błąd_Względny_Rozwiązania', ...
    'Współczynnik_Stabilności', ...
    'Współczynnik_Poprawności'}));
