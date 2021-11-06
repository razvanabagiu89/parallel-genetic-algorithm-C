# Parallel Genetic Algorithm in C

## About

https://curs.upb.ro/2021/pluginfile.php/367342/mod_resource/content/6/Tema%201%20-%20Enunt.pdf

Compute a parallelisation of a genetic algorithm used as an example of Darwin's theory based on the fittest -> the more likely to survive

## How to run the tests

- run test script local or in Docker dir

## Structure

### For solving the tasks I used the following files:

- tema1_par.c 

The whole computation part is done in this file with the using of the <pthread.h> library

- sack_object.h & individual.h & genetic_algorithm.h

Same as skel files with slight modifications


## Flow and parallelisation

- main

Contains parsing of the input and allocation of all the pointers used. By using the
<pthread.h> library here the threads are created and used in the 'run_genetic_algorithm' function, then all threads are joined and the memory is freed. Return is checked on every function call to ensure reliability of the program.

- run_genetic_algorithm

All functions from the skel are parallelised alongside with the sorting algorithm. First of all, the generations are initialised and their chromosomes are allocated, then each thread computes fitness for its part of the array. I used barrier to ensure the synchronization of the threads, and to prevent the race conditions. In some parts of the code, only thread[0] (Work-Crew design) does the work so the other threads don't take part and facilitate the race condition. Then, each thread sorts its part of the array by using qsort, and then thread[0] merges all the parts into one within a while() loop. All the parts are being sorted at the same time with qsort which ensures scalability and efficiency. Then, computation of the elites, bit string mutation one and two and crossover are done using paralleled arrays and skel functions. The same fitness computation and sorting applies after the generations loop is finished. All the resources are freed on the final part of the program.



## Project structure including tests
```bash
├── Makefile
├── genetic_algorithm.h
├── individual.h
├── inputs
│   ├── in0
│   ├── in1
│   ├── in2
│   ├── in3
│   └── in4
├── sack_object.h
├── tema1_par
└── tema1_par.c


```
