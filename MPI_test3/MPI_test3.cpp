//подключаем три библиотеки
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <ctime>
#include <fstream>
#include <vector>
#include <string>
#define N 10000

using namespace std;

const int arr_size = 10000;

void num_file()
{
	ofstream out;

	out.open("numbers1.txt");
	for (int i = 0; i < arr_size; i++)
	{
		for (int j = 0; j < arr_size; j++)
		{
			out << 1 + rand() % 10 << " ";
		}
		out << endl;
	}
	out.close();

	out.open("numbers2.txt");
	for (int i = 0; i < arr_size; i++)
	{
		for (int j = 0; j < arr_size; j++)
		{
			out << 1 + rand() % 10 << " ";
		}
		out << endl;
	}
	out.close();
}

int main(int argc, char** argv)
{
	int rank, size;
	MPI_Status status;
	double sum = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size); //кол-во запущен процессов (size)
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//-----------------------------------------------------------//
	//srand(time(0));
	//if (rank == 0) num_file();

	double* arr1 = new double[arr_size*arr_size];
	double* arr2 = new double[arr_size*arr_size];

	//unsigned int start_time = clock();
	double t1 = MPI_Wtime();

	// Считаем матрицу из файла
	if (rank == 0)
	{
		ifstream in("numbers1.txt");
		for (int i = 0; i < arr_size; i++)
			for (int j = 0; j < arr_size; j++)
				in >> arr1[i*arr_size + j];
		in.close();
		
		ifstream in2;
		in2.open("numbers2.txt");
		for (int i = 0; i < arr_size; i++)
			for (int j = 0; j < arr_size; j++)
				in2 >> arr2[i*arr_size+j];
		in2.close();
	}
	
	MPI_Bcast(&arr1[0], arr_size * arr_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&arr2[0], arr_size * arr_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	int k = arr_size / size;
	int start = rank * k; // если 1 процесс то стартовая позиция = 25
	int stop = (rank + 1) * k;

	if (rank == size - 1) stop = arr_size;
	// Евклидова норма матрицы
	double S = 0;
	for (int i = start; i < stop; ++i)
	{
		for (int j = 0; j < arr_size; ++j)
		{
			S += arr1[i*arr_size + j] * arr1[i*arr_size + j];
		}
	}
	
	double total = 0;

	MPI_Reduce(&S, &total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&S, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (rank == 0)
	{
		total = sqrt(total);
		printf("E = %f\n", total);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	double* arr3 = new double[arr_size*arr_size];

	for (int i = start; i < stop; ++i)
	{
		for (int j = 0; j < arr_size; ++j)
		{
			arr3[i*arr_size + j] = total * arr1[i*arr_size + j] + arr2[i*arr_size + j] - arr1[i*arr_size + j];
		}
	}

	if (rank == 0)
	{
		ofstream out("out.txt");
		for (int i = 0; i < arr_size; i++)
		{
			for (int j = 0; j < arr_size; j++)
			{
				out << arr3[i*arr_size + j] << " ";
			}
			out << endl;
		}
	}

	double t2 = MPI_Wtime() - t1;

	if (rank == 0) printf("time = %f\n", t2);

	MPI_Finalize();
	
	return(0);
}
