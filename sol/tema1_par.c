#include "genetic_algorithm.h"

// global variables for ease of access

// sack of objects
sack_object *objects;
// number of objects that can be put in the sack
int object_count;
// max capacity of the sack
int sack_capacity;
// number of generations to be computed
int generations_count;
// number of threads
int P;
// each thread's id
int id;

// current computed generation
individual *current_generation;
// next generation to be computed
individual *next_generation;
// sorted generation after the sorting algorithm is performed
individual *sorted_generation;
// aux used for swapping the generations
individual *aux;
// array used for merge algorithm to remember last end position of the thread -> [0] reference
int *thread_end;

// barrier used for synchronization of the threads
pthread_barrier_t barrier;

// used for merging all the sorted parts done by the threads
void merge(individual *generation, int left, int mid, int right)
{
	int size1 = mid - left + 1;
	int size2 = right - mid;

	// temp arrays for storing each sorted part
	individual first[size1], second[size2];

	for (int i = 0; i < size1; i++)
	{
		first[i] = generation[left + i];
	}
	for (int j = 0; j < size2; j++)
	{
		second[j] = generation[mid + 1 + j];
	}

	int i = 0;
	int j = 0;
	int k = left;

	while (i < size1 && j < size2)
	{
		// first sorting criteria: descending fitness
		if (first[i].fitness > second[j].fitness)
		{
			memcpy(generation + k, first + i, sizeof(individual));
			i++;
		}
		else if (first[i].fitness < second[j].fitness)
		{
			memcpy(generation + k, second + j, sizeof(individual));
			j++;
		}
		// if fitness is equal then second sorting criteria: ascending number of objects
		else
		{
			int first_count = 0, second_count = 0;

			for (int k = 0; k < first[i].chromosome_length && k < second[j].chromosome_length; ++k)
			{
				first_count += first[i].chromosomes[k];
				second_count += second[j].chromosomes[k];
			}

			if (first_count < second_count)
			{
				memcpy(generation + k, first + i, sizeof(individual));
				i++;
			}
			else if (first_count > second_count)
			{
				memcpy(generation + k, second + j, sizeof(individual));
				j++;
			}
			// if number of objects is equal then third sorting criteria: descending index
			else
			{
				if (first[i].index > second[j].index)
				{
					memcpy(generation + k, first + i, sizeof(individual));
					i++;
				}
				else if (first[i].index <= second[j].index)
				{
					memcpy(generation + k, second + j, sizeof(individual));
					j++;
				}
			}
		}
		k++;
	}

	// if one array is finished, copy the rest of elements
	while (i < size1)
	{
		memcpy(generation + k, first + i, sizeof(individual));
		i++;
		k++;
	}

	while (j < size2)
	{
		memcpy(generation + k, second + j, sizeof(individual));
		j++;
		k++;
	}
}

// skel function
void print_best_fitness(const individual *generation)
{
	printf("%d\n", generation[0].fitness);
}

// used for debugging purposes
void print_generation(const individual *generation, int limit)
{
	for (int i = 0; i < limit; ++i)
	{
		for (int j = 0; j < generation[i].chromosome_length; ++j)
		{
			printf("%d ", generation[i].chromosomes[j]);
		}
		printf("\n%d - %d\n", i, generation[i].fitness);
	}
}

// skel function
void free_generation(individual *generation)
{
	int i;

	for (i = 0; i < generation->chromosome_length; ++i)
	{
		free(generation[i].chromosomes);
		generation[i].chromosomes = NULL;
		generation[i].fitness = 0;
	}
}

// skel function
void mutate_bit_string_2(const individual *ind, int generation_index)
{
	int step = 1 + generation_index % (ind->chromosome_length - 2);

	// mutate all chromosomes by a given step
	for (int i = 0; i < ind->chromosome_length; i += step)
	{
		ind->chromosomes[i] = 1 - ind->chromosomes[i];
	}
}

// skel function
void crossover(individual *parent1, individual *child1, int generation_index)
{
	individual *parent2 = parent1 + 1;
	individual *child2 = child1 + 1;
	int count = 1 + generation_index % parent1->chromosome_length;

	memcpy(child1->chromosomes, parent1->chromosomes, count * sizeof(int));
	memcpy(child1->chromosomes + count, parent2->chromosomes + count, (parent1->chromosome_length - count) * sizeof(int));

	memcpy(child2->chromosomes, parent2->chromosomes, count * sizeof(int));
	memcpy(child2->chromosomes + count, parent1->chromosomes + count, (parent1->chromosome_length - count) * sizeof(int));
}

// skel function
void mutate_bit_string_1(const individual *ind, int generation_index)
{
	int i, mutation_size;
	int step = 1 + generation_index % (ind->chromosome_length - 2);

	if (ind->index % 2 == 0)
	{
		// for even-indexed individuals, mutate the first 40% chromosomes by a given step
		mutation_size = ind->chromosome_length * 4 / 10;
		for (i = 0; i < mutation_size; i += step)
		{
			ind->chromosomes[i] = 1 - ind->chromosomes[i];
		}
	}
	else
	{
		// for even-indexed individuals, mutate the last 80% chromosomes by a given step
		mutation_size = ind->chromosome_length * 8 / 10;
		for (i = ind->chromosome_length - mutation_size; i < ind->chromosome_length; i += step)
		{
			ind->chromosomes[i] = 1 - ind->chromosomes[i];
		}
	}
}

// skel function
void copy_individual(const individual *from, const individual *to)
{
	memcpy(to->chromosomes, from->chromosomes, from->chromosome_length * sizeof(int));
}

// skel function
int cmpfunc(const void *a, const void *b)
{
	int i;
	individual *first = (individual *)a;
	individual *second = (individual *)b;

	int res = second->fitness - first->fitness; // decreasing by fitness
	if (res == 0)
	{
		int first_count = 0, second_count = 0;

		for (i = 0; i < first->chromosome_length && i < second->chromosome_length; ++i)
		{
			first_count += first->chromosomes[i];
			second_count += second->chromosomes[i];
		}

		res = first_count - second_count; // increasing by number of objects in the sack
		if (res == 0)
		{
			return second->index - first->index;
		}
	}

	return res;
}

// parallelised algorithm function
void *run_genetic_algorithm(void *arg)
{

	int thread_id = *(int *)arg;
	int start = thread_id * (int)ceil(1.0f * object_count / P);
	int end = MIN(object_count, (thread_id + 1) * (int)ceil(1.0f * object_count / P));

	// paralleled initial generation (composed of object_count individuals with a single item in the sack)
	for (int i = start; i < end; ++i)
	{
		current_generation[i].fitness = 0;
		current_generation[i].chromosomes = (int *)calloc(object_count, sizeof(int));
		current_generation[i].chromosomes[i] = 1;
		current_generation[i].index = i;
		current_generation[i].chromosome_length = object_count;

		next_generation[i].fitness = 0;
		next_generation[i].chromosomes = (int *)calloc(object_count, sizeof(int));
		next_generation[i].index = i;
		next_generation[i].chromosome_length = object_count;

		sorted_generation[i].chromosomes = (int *)calloc(object_count, sizeof(int));
	}

	int count, cursor;
	start = thread_id * (int)ceil(1.0f * object_count / P);
	end = MIN(object_count, (thread_id + 1) * (int)ceil(1.0f * object_count / P));

	for (int k = 0; k < generations_count; ++k)
	{
		cursor = 0;
		int weight, profit;

		// paralleled fitness computation
		for (int i = start; i < end; ++i)
		{
			weight = 0;
			profit = 0;

			for (int j = 0; j < current_generation[i].chromosome_length; ++j)
			{
				if (current_generation[i].chromosomes[j])
				{
					weight += objects[j].weight;
					profit += objects[j].profit;
				}
			}
			current_generation[i].fitness = (weight <= sack_capacity) ? profit : 0;
		}
		pthread_barrier_wait(&barrier);

		// [0] paralleled sort algorithm - qsort + merge sort
		start = thread_id * (int)ceil(1.0f * object_count / P);
		end = MIN(object_count, (thread_id + 1) * (int)ceil(1.0f * object_count / P));

		if (thread_id == (P - 1) && thread_id != 0 && object_count % P != 0)
		{
			end = object_count;
		}
		thread_end[thread_id] = end;
		int nr_elems = end - start;

		// each thread sorts its own part of the array
		qsort(current_generation + start, nr_elems, sizeof(individual), cmpfunc);
		pthread_barrier_wait(&barrier);

		int nr_merges = P - 1;
		if (thread_id == 0 && P > 1)
		{
			// initialization of indices
			int left = 0;
			int id = 1;
			end = thread_end[id] - 1;
			int mid = end / 2;

			while (nr_merges > 0)
			{
				merge(current_generation, left, mid, end);
				nr_merges--;
				thread_end[id] = end;
				id++;
				end = thread_end[id] - 1;
				mid = thread_end[id - 1];
			}
		}
		pthread_barrier_wait(&barrier);

		// keep first 30% children (elite children selection)
		count = object_count * 3 / 10;
		start = thread_id * (int)ceil(1.0f * count / P);
		end = MIN(count, (thread_id + 1) * (int)ceil(1.0f * count / P));

		// paralleled copy
		for (int i = start; i < end; ++i)
		{
			copy_individual(current_generation + i, next_generation + i);
		}

		if (thread_id == 0)
		{
			cursor += count;
		}
		pthread_barrier_wait(&barrier);

		// mutate first 20% children with the first version of bit string mutation
		count = object_count * 2 / 10;
		start = thread_id * (int)ceil(1.0f * count / P);
		end = MIN(count, (thread_id + 1) * (int)ceil(1.0f * count / P));

		// paralleled copy and bit string mutate 1
		for (int i = start; i < end; ++i)
		{
			copy_individual(current_generation + i, next_generation + cursor + i);
			mutate_bit_string_1(next_generation + cursor + i, k);
		}

		if (thread_id == 0)
		{
			cursor += count;
		}
		pthread_barrier_wait(&barrier);

		// mutate next 20% children with the second version of bit string mutation
		count = object_count * 2 / 10;
		start = thread_id * (int)ceil(1.0f * count / P);
		end = MIN(count, (thread_id + 1) * (int)ceil(1.0f * count / P));

		// paralleled copy and bit string mutate 2
		for (int i = start; i < end; ++i)
		{
			copy_individual(current_generation + i + count, next_generation + cursor + i);
			mutate_bit_string_2(next_generation + cursor + i, k);
		}

		if (thread_id == 0)
		{
			cursor += count;
		}
		pthread_barrier_wait(&barrier);

		// crossover first 30% parents with one-point crossover
		// (if there is an odd number of parents, the last one is kept as such)
		count = object_count * 3 / 10;
		pthread_barrier_wait(&barrier);

		// paralleled crossover
		if (thread_id == 0)
		{
			if (count % 2 == 1)
			{
				copy_individual(current_generation + object_count - 1, next_generation + cursor + count - 1);
				count--;
			}
		}
		pthread_barrier_wait(&barrier);

		start = thread_id * (int)ceil(1.0f * count / P);
		end = MIN(count, (thread_id + 1) * (int)ceil(1.0f * count / P));

		for (int i = start; i < end; i += 2)
		{
			crossover(current_generation + i, next_generation + cursor + i, k);
		}
		pthread_barrier_wait(&barrier);

		// swap generations
		if (thread_id == 0)
		{
			aux = current_generation;
			current_generation = next_generation;
			next_generation = aux;
		}

		start = thread_id * (int)ceil(1.0f * object_count / P);
		end = MIN(object_count, (thread_id + 1) * (int)ceil(1.0f * object_count / P));

		// reset indices
		for (int i = start; i < end; ++i)
		{
			current_generation[i].index = i;
		}
		pthread_barrier_wait(&barrier);

		if (thread_id == 0)
		{
			if (k % 5 == 0)
			{
				print_best_fitness(current_generation);
			}
		}
	}

	// same paralleled algorithm done for fitness + sort
	// compute fitness and sort by it
	start = thread_id * (int)ceil(1.0f * object_count / P);
	end = MIN(object_count, (thread_id + 1) * (int)ceil(1.0f * object_count / P));

	int weight = 0;
	int profit = 0;

	for (int i = start; i < end; ++i)
	{
		weight = 0;
		profit = 0;

		for (int j = 0; j < current_generation[i].chromosome_length; ++j)
		{
			if (current_generation[i].chromosomes[j])
			{
				weight += objects[j].weight;
				profit += objects[j].profit;
			}
		}

		current_generation[i].fitness = (weight <= sack_capacity) ? profit : 0;
	}
	pthread_barrier_wait(&barrier);

	start = thread_id * (int)ceil(1.0f * object_count / P);
	end = MIN(object_count, (thread_id + 1) * (int)ceil(1.0f * object_count / P));

	if (thread_id == (P - 1) && thread_id != 0 && object_count % P != 0)
	{
		end = object_count;
	}

	thread_end[thread_id] = end;
	int nr_elems = end - start;

	qsort(current_generation + start, nr_elems, sizeof(individual), cmpfunc);
	pthread_barrier_wait(&barrier);

	int nr_merges = P - 1;
	if (thread_id == 0 && P > 1)
	{
		// initialization of indices
		int left = 0;
		int id = 1;
		end = thread_end[id] - 1;
		int mid = end / 2;

		while (nr_merges > 0)
		{
			merge(current_generation, left, mid, end);
			nr_merges--;
			thread_end[id] = end;
			id++;
			end = thread_end[id] - 1;
			mid = thread_end[id - 1];
		}
	}

	pthread_barrier_wait(&barrier);

	if (thread_id == 0)
	{
		print_best_fitness(current_generation);

		// free resources for old generation
		free_generation(current_generation);
		free_generation(next_generation);

		// free resources
		free(current_generation);
		free(next_generation);
	}

	return NULL;
}

int main(int argc, char *argv[])
{
	if (argc < 4)
	{
		fprintf(stderr, "Usage:\n\t./tema1_par in_file generations_count P\n");
		return 0;
	}

	object_count = 0;
	sack_capacity = 0;
	generations_count = 0;
	P = 0;

	// input parse
	FILE *fp;
	fp = fopen(argv[1], "r");
	if (fp == NULL)
	{
		return 0;
	}

	if (fscanf(fp, "%d %d", &object_count, &sack_capacity) < 2)
	{
		fclose(fp);
		return 0;
	}

	if (object_count % 10)
	{
		fclose(fp);
		return 0;
	}

	objects = (sack_object *)calloc(object_count, sizeof(sack_object));

	for (int i = 0; i < object_count; ++i)
	{
		if (fscanf(fp, "%d %d", &objects[i].profit, &objects[i].weight) < 2)
		{
			free(objects);
			fclose(fp);
			return 0;
		}
	}

	fclose(fp);

	generations_count = atoi(argv[2]);

	if (generations_count == 0)
	{
		free(objects);

		return 0;
	}

	P = atoi(argv[3]);

	int thread_id[P];
	pthread_t tid[P];
	int r;

	pthread_barrier_init(&barrier, NULL, P);

	current_generation = (individual *)calloc(object_count, sizeof(individual));
	next_generation = (individual *)calloc(object_count, sizeof(individual));
	sorted_generation = (individual *)calloc(object_count, sizeof(individual));
	thread_end = (int *)calloc(P, sizeof(int));

	if (current_generation == NULL || next_generation == NULL || sorted_generation == NULL || thread_end == NULL)
	{
		printf("Error on malloc!\n");
		fflush(stdout);
		exit(1);
	}

	for (int i = 0; i < P; i++)
	{
		thread_id[i] = i;
		r = pthread_create(&tid[i], NULL, run_genetic_algorithm, &thread_id[i]);

		if (r)
		{
			printf("Error on thread creation, thread id: %d\n", i);
			exit(-1);
		}
	}

	for (int i = 0; i < P; i++)
	{
		r = pthread_join(tid[i], NULL);

		if (r)
		{
			printf("Error on waiting thread: %d\n", i);
			exit(-1);
		}
	}
	pthread_barrier_destroy(&barrier);

	free(objects);
	return 0;
}
