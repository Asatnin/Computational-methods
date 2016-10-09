function Lab1_Hungarian_alg()
clc
c = read_matrix();
init_c = c;
n = size(c, 1);
c = reduce_matrix_columns(c);
c = reduce_matrix_rows(c);
[z, k_zeros] = init_independent_sys_zeros(c);
marked_rows = [];
need_mark_columns = true;

while k_zeros ~= n
    if need_mark_columns
        marked_columns = mark_columns(z);
    end
    [z, has_prime_zero, col, row] = mark_prime_zero(c, z, marked_columns, ...
        marked_rows);
    while ~has_prime_zero
        c = prepare_prime_zero(c, marked_columns, marked_rows);
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
        z = build_L_chain(z, col, row);
        marked_rows = [];
        marked_columns = [];
        k_zeros = k_zeros + 1;
    end
end
print_answer(init_c, z);
%[z, has_prime_zero, col, row] = mark_prime_zero(c, z, marked_columns, marked_rows)
end

function c = read_matrix()
%c = [9 11 3 6 6; 10 9 11 5 6; 8 10 5 6 4; 6 8 10 4 9; 11 10 9 8 7];
c = importdata('input.txt')
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

function [res, has_prime_zero, col, row] = mark_prime_zero(c, z, ...
    marked_columns, marked_rows)
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

function res = prepare_prime_zero(c, marked_columns, marked_rows)
h = 0;
has_min = false;
for j = 1:size(c, 2)
    if ~ismember(j, marked_columns)
        for i = 1:size(c, 1)
            if (~ismember(i, marked_rows)) && (~hasmin || (c(i,j) < h))
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

res = c;
end

function res = build_L_chain(z, col, row)
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
    else
        col = next_horizontal(z, row);
        if col == -1
            break
        end
        res(row,col) = 1;
    end
    step = step + 1;
end
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
x
disp('Значение целевой функции');
f
end