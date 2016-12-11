function Lab1_Hungarian_alg()
clc
%filename = input('Введите название файла с исходными данными: ', 's');
filename = 'input.txt';
problem_type = input(['Выберите тип задачи (1 - минимизация, ' ...
    '2 - максимизация): ']);
debug = input('Выберите режим работы (1 - итоговый, 2 - отладочный): ');
c = read_matrix(filename);
init_c = c;
n = size(c, 1);
if problem_type == 2 % максимизация
    c = -c + max(max(c));
end
c = reduce_matrix_columns(c);
c = reduce_matrix_rows(c);
[z, k_zeros] = init_independent_sys_zeros(c);
marked_rows = [];
need_mark_columns = true;

% номер итерации (для debug-режима)
it = 0;
if debug == 2
    %print_init_c(init_c);
    disp(strcat('Итерация #', int2str(it)));
    print_after_init_c(c);
    %disp('Исходная СНН:');
    print_cur_iteration(c, z);
end

% основной этап венгерского метода
while k_zeros ~= n
    flag = false;
    h = -1;
    
    if need_mark_columns
        marked_columns = mark_columns(z);
    end
    [z, has_prime_zero, col, row] = mark_prime_zero(c, z, marked_columns, ...
        marked_rows);
    while ~has_prime_zero
        [c, h] = prepare_prime_zero(c, marked_columns, marked_rows);
        flag = true;
        [z, has_prime_zero, col, row] = mark_prime_zero(c, z, ...
            marked_columns, marked_rows);
    end
    
    need_mark_columns = true;
    for j = 1:size(z, 2)
        if (j ~= col) && (z(row, j) == 1)
            marked_columns(marked_columns == j) = [];
            marked_rows(end + 1) = row;
            need_mark_columns = false;
            break;
        end
    end
    
    if need_mark_columns
        prev_z = z;
        [z, rows, cols] = build_L_chain(z, col, row);
        marked_rows = [];
        marked_columns = [];
        k_zeros = k_zeros + 1;
        
        if debug == 2
            it = it + 1;
            disp(strcat('Итерация #', int2str(it)));
            if flag
                print_new_c(c, h);
            end
            print_L_chain(rows, cols);
            print_cur_iteration2(c, prev_z);
            print_cur_iteration(c, z);
        end
    end
end
print_answer(init_c, z);
end

% чтение исходных данных из файла
function c = read_matrix(filename)
c = importdata(filename);
end

function res = reduce_matrix_columns(c)
for j = 1:size(c, 2)
    c(:,j) = c(:,j) - min(c(:,j));
end
res = c;
end

function res = reduce_matrix_rows(c)
for i = 1:size(c, 1)
    c(i,:) = c(i,:) - min(c(i,:));
end
res = c;
end

% построение первоначальной системы независимых нулей
function [z, k_zeros, marked_columns] = init_independent_sys_zeros(c)
marked_columns = [];
k_zeros = 0;
z = zeros(size(c));
for j = 1:size(c, 2)
    for i = 1:size(c, 1)
        break_loop = false;
        if c(i,j) == 0
            independent = true;            
            for k = 1:size(c, 2)
                if (j ~= k) && (z(i,k) == 1)
                    independent = false;
                    break
                end
            end
            
            if independent
                z(i,j) = 1;
                k_zeros = k_zeros + 1;
                break_loop = true;
            end
        end
        
        if break_loop
            break
        end
    end
end
end

function res = mark_columns(z)
res = [];
for j = 1:size(z, 2)
    for i = 1:size(z, 1)
        if z(i,j) == 1
            res(end + 1) = j;
            break;
        end
    end
end
end

% поиск и пометка нуля среди невыделенных элементов
function [res, has_prime_zero, col, row] = mark_prime_zero(c, z, ...
    marked_columns, marked_rows)
col = -1;
row = -1;
for j = 1:size(c, 2)
    if ~ismember(j, marked_columns)
        break_loop = false;
        for i = 1:size(c, 1)
            if (~ismember(i, marked_rows)) && (c(i,j) == 0)
                col = j;
                row = i;
                z(i,j) = 2;
                break_loop = true;
                break
            end
        end
        if break_loop
            break
        end
    end
end

res = z;
has_prime_zero = break_loop;
end

function [res, h_ans] = prepare_prime_zero(c, marked_columns, marked_rows)
h = 0;
has_min = false;
for j = 1:size(c, 2)
    if ~ismember(j, marked_columns)
        for i = 1:size(c, 1)
            if (~ismember(i, marked_rows)) && (~has_min || (c(i,j) < h))
                has_min = true;
                h = c(i,j);
            end
        end
    end
end

for i = 1:size(c, 1)
    if ~ismember(i, marked_rows)
c(i,:) = c(i,:) - h;
    end
end

for j = 1:size(c, 2)
    if ismember(j, marked_columns)
        c(:,j) = c(:,j) + h;
    end
end

h_ans = h;
res = c;
end

% построение L-цепочки
function [res, rows_ans, cols_ans] = build_L_chain(z, col, row)
rows = [];
cols = [];
rows(end + 1) = row;
cols(end + 1) = col;
res = z;
res(row,col) = 1;
step = 0;
while true
    if mod(step, 2) == 0
        row = next_vertical(z, col);
        if row == -1
            break
        end
        res(row,col) = 0;
        rows(end + 1) = row;
        cols(end + 1) = col;
    else
        col = next_horizontal(z, row);
        if col == -1
            break
        end
        res(row,col) = 1;
        rows(end + 1) = row;
        cols(end + 1) = col;
    end
    step = step + 1;
end
res(res == 2) = 0;
rows_ans = rows;
cols_ans = cols;
end

function row = next_vertical(z, col)
row = -1;
for i = 1:size(z, 1)
    if z(i,col) == 1
        row = i;
        break
    end
end
end

function col = next_horizontal(z, row)
col = -1;
for j = 1:size(z, 2)
    if z(row, j) == 2
        col = j;
        break
    end
end
end

function print_answer(c, z)
disp('Искомая матрица назначений:');
f = 0;
x = zeros(size(z));
for i = 1:size(c, 1)
    for j = 1:size(c, 2)
        if z(i,j) == 1
            x(i,j) = 1;
            f = f + c(i,j);
        end
    end
end
disp(x);
disp('Значение целевой функции');
disp(f);
end

% вывод отладочной информации - СНН на каждом шаге
function print_cur_iteration(c, z)
disp('Текущая СНН:');
for i = 1:size(c, 1)
    fprintf('     ');
    for j = 1:size(c, 2)
        if z(i,j) == 0
            fprintf('%d     ', c(i,j));
        elseif z(i,j) == 1
            fprintf('%d*    ', c(i,j));
        else
            fprintf('%d''    ', c(i,j));
        end
    end
    fprintf('\n');
end
fprintf('\n\n');
end

function print_cur_iteration2(c, z)
for i = 1:size(c, 1)
    fprintf('     ');
    for j = 1:size(c, 2)
        if z(i,j) == 0
            fprintf('%d     ', c(i,j));
        elseif z(i,j) == 1
            fprintf('%d*    ', c(i,j));
        else
            fprintf('%d''    ', c(i,j));
        end
    end
    fprintf('\n');
end
fprintf('\n');
end

function print_new_c(c, h)
fprintf('h = %d\n', h);
disp('Преобразованная матрица стоимостей: ');
for i = 1:size(c, 1)
    fprintf('     ');
    for j = 1:size(c, 2)
        fprintf('%d     ', c(i,j));
    end
    fprintf('\n');
end
end

function print_init_c(c)
disp('Исходная матрица стоимостей: ');
for i = 1:size(c, 1)
    fprintf('     ');
    for j = 1:size(c, 2)
        fprintf('%d     ', c(i,j));
    end
    fprintf('\n');
end
end

function print_after_init_c(c)
disp('Преобразованная матрица стоимостей: ');
for i = 1:size(c, 1)
    fprintf('     ');
    for j = 1:size(c, 2)
        fprintf('%d     ', c(i,j));
    end
    fprintf('\n');
end
end

function print_L_chain(rows, cols)
fprintf('Построенная L-цепочка состоит из %d элементов:\n', size(rows, 2));
for i = 1:(size(rows, 2) - 1)
    fprintf('(%d, %d) --> ', rows(i), cols(i));
end
fprintf('(%d, %d)\n', rows(end), cols(end));
end
