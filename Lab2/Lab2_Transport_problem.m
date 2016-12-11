function Lab2_Transport_problem()
clc
[s, d, c] = read_from_file();
c_init = c;
%problem_type = input(['Выберите тип задачи (1 - минимизация, ' ...
%    '2 - максимизация): ']);
%if problem_type == 2 % максимизация
%    c = -c + max(max(c));
%end
[x, bas] = north_west_corner(s, d);
print_init_solution(x);
it = 0;
while true
    [u, v] = make_and_solve_equations_for_basis_cells(bas, c);
    [row, col] = find_cell_to_include_in_basis(u, v, bas, c);
    if (row == -1)
        break;
    end
    bas(row, col) = 1;
    [cycle_rows, cycle_cols] = find_cycle(row, col, bas);
    it = it + 1;
    print_cur_iteration(x, cycle_rows, cycle_cols, it)
    [x, bas] = reallocate(x, bas, cycle_rows, cycle_cols);
end
print_final_solution(x, c_init);
end


function print_init_solution(x)
disp('Начальное БДР: ');
for i = 1:size(x, 1)
    for j = 1:size(x, 2)
        fprintf('%6d', x(i,j));
    end
    fprintf('\n');
end
end


function print_cur_iteration(x, rows, cols, it)
fprintf('\n**************************\n')
fprintf('Итерация #%d.\n\nТекущая транспортная таблица:\n', it);
for i = 1:size(x, 1)
    for j = 1:size(x, 2)
        fprintf('%6d', x(i,j));
    end
    fprintf('\n');
end
fprintf('\nТранспортный цикл состоит из %d элементов:\n', size(rows, 1));
for i = 1:(size(rows, 1) - 1)
    fprintf('(%d, %d) --> ', rows(i), cols(i));
end
fprintf('(%d, %d)\n', rows(end), cols(end));
end


function print_final_solution(x, c)
fprintf('\n**************************\n')
fprintf('Итоговая транспортная таблица:\n');
sum = 0;
for i = 1:size(x, 1)
    for j = 1:size(x, 2)
        sum = sum + x(i, j) * c(i, j);
        fprintf('%6d', x(i,j));
    end
    fprintf('\n');
end
fprintf('\nЗначение целевой функции: %d\n', sum);
end


function [x_new, bas_new] = reallocate(x, bas, cycle_rows, cycle_cols)
n = size(cycle_rows, 1);
w = x(cycle_rows(2), cycle_cols(2));
row = -1;
col = -1;
for i = 2:2:n
    if (x(cycle_rows(i), cycle_cols(i)) < w)
        w = x(cycle_rows(i), cycle_cols(i));
    end
end
for i = 1:n
    if (mod(i, 2) == 1)
        x(cycle_rows(i), cycle_cols(i)) = x(cycle_rows(i), cycle_cols(i)) + w;
    else
        x(cycle_rows(i), cycle_cols(i)) = x(cycle_rows(i), cycle_cols(i)) - w;
        if (x(cycle_rows(i), cycle_cols(i)) == 0)
            row = cycle_rows(i);
            col = cycle_cols(i);
        end
    end
end
bas(row, col) = 0;
x_new = x;
bas_new = bas;
end



function [cycle_rows, cycle_cols] = find_cycle(init_row, init_col, bas)
m = size(bas, 1);
n = size(bas, 2);
cycle_rows = []; cycle_rows(end + 1) = init_row;
cycle_cols = []; cycle_cols(end + 1) = init_col;
[cycle_rows, cycle_cols, good] = look_in_row(cycle_rows, cycle_cols, ...
    init_row, init_col, init_row, init_col, bas); 
end


function [cycle_rows, cycle_cols, good] = look_in_column(prev_cycle_rows, ...
    prev_cycle_cols, row, col, init_row, init_col, bas)
m = size(bas, 1);
for i = 1:m
    if (i ~= row && bas(i, col) == 1)
        [new_rows, new_cols, good_row] = look_in_row(prev_cycle_rows, ...
            prev_cycle_cols, i, col, init_row, init_col, bas);
        if (good_row)
            cycle_rows = [new_rows; i];
            cycle_cols = [new_cols; col];
            good = true;
            return;
        end
    end
end
cycle_rows = prev_cycle_rows;
cycle_cols = prev_cycle_cols;
good = false;
end


function [cycle_rows, cycle_cols, good] = look_in_row(prev_cycle_rows, ...
    prev_cycle_cols, row, col, init_row, init_col, bas)
n = size(bas, 2);
for j = 1:n
    if (j ~= col && bas(row, j) == 1)
        if (j == init_col)
            cycle_rows = [prev_cycle_rows; row];
            cycle_cols = [prev_cycle_cols; j];
            good = true;
            return;
        end
        
        [new_rows, new_cols, good_col] = look_in_column(prev_cycle_rows, ...
            prev_cycle_cols, row, j, init_row, init_col, bas);
        if (good_col)
            cycle_rows = [new_rows; row];
            cycle_cols = [new_cols; j];
            good = true;
            return;
        end
    end
end
cycle_rows = prev_cycle_rows;
cycle_cols = prev_cycle_cols;
good = false;
end


function [row, col] = find_cell_to_include_in_basis(u, v, bas, c)
m = size(c, 1);
n = size(c, 2);
row = -1;
col = -1;
d_min = 0;
for i = 1:m
    for j = 1:n
        if (bas(i, j) == 0 && c(i, j) - u(i) - v(j) < d_min)
            d_min = c(i, j) - u(i) - v(j);
            row = i;
            col = j;
        end
    end
end
end


function [u, v] = make_and_solve_equations_for_basis_cells(bas, c)
m = size(c, 1);
n = size(c, 2);
used_u = zeros(m, 1);
used_v = zeros(n, 1);
u = zeros(m, 1);
v = zeros(n, 1);
used_u(1) = 1;
num_of_unknowns = 1;
while (num_of_unknowns ~= m + n)
    for i = 1:m
        for j = 1:n
            if (bas(i, j) == 1)
                if (used_u(i) == 0 && used_v(j) == 1)
                    u(i) = c(i, j) - v(j);
                    num_of_unknowns = num_of_unknowns + 1;
                    used_u(i) = 1;
                end
                if (used_u(i) == 1 && used_v(j) == 0)
                    v(j) = c(i, j) - u(i);
                    num_of_unknowns = num_of_unknowns + 1;
                    used_v(j) = 1;
                end
            end
        end
    end
end
end


function [x, bas] = north_west_corner(s, d)
m = size(s, 2);
n = size(d, 2);
x = zeros(m, n);
bas = zeros(m, n);
i = 1;
j = 1;
while ((i <= m) && (j <= n))
    if (d(j) <= s(i))  % для вырожденного случая вычеркнем столбец
        bas(i, j) = 1;
        x(i, j) = d(j);
        s(i) = s(i) - d(j);
        j = j + 1;
    else
        bas(i, j) = 1;
        x(i, j) = s(i);
        d(j) = d(j) - s(i);
        i = i + 1;
    end
end
end


function [s, d, c] = read_from_file()
filename = 'input.txt';
fid = fopen(filename);
s_line = fgets(fid);
s = str2num(s_line);
d_line = fgets(fid);
d = str2num(d_line);
c = [];
for i = 1:size(s, 2)
    tmp_line = fgets(fid);
    c = [c; str2num(tmp_line)];
end
fclose(fid);
end