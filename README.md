# Параллелизация с помощью директивы for

## Описание фаилов
- var40.c - файл с исходным кодом
- 

## Полезная информация
Авторизация на Polus происходит через установление SSH-соединения

1. Для компиляции фаила на вычислительном кластере Polus: 
gcc -std=c99 -fopenmp var40_variations.c
или 
xlc -qsmp=omp -qarch=pwr8 for.c 

2. Для запуска фаила 
mpisubmit.pl --stdout=results.txt ./a.out
