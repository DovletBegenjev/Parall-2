//подключаем три библиотеки
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <ctime>
#include <fstream>
#include <vector>
#include <string>
#define N 10

using namespace std;

const int arr_size = 1000;

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
		ifstream in("numbers100.txt");
		for (int i = 0; i < arr_size; i++)
			for (int j = 0; j < arr_size; j++)
				in >> arr1[i*arr_size + j];
		in.close();
		
		ifstream in2;
		in2.open("numbers200.txt");
		for (int i = 0; i < arr_size; i++)
			for (int j = 0; j < arr_size; j++)
				in2 >> arr2[i*arr_size+j];
		in2.close();
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0)
	{
		for (int i = 1; i < size; i++)
		{
			MPI_Send(&arr1[0], arr_size * arr_size, MPI_DOUBLE, i, 77, MPI_COMM_WORLD);
			MPI_Send(&arr2[0], arr_size * arr_size, MPI_DOUBLE, i, 55, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Recv(&arr1[0], arr_size * arr_size, MPI_DOUBLE, 0, 77, MPI_COMM_WORLD, &status);
		MPI_Recv(&arr2[0], arr_size * arr_size, MPI_DOUBLE, 0, 55, MPI_COMM_WORLD, &status);
	}
	//MPI_Barrier(MPI_COMM_WORLD);
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
	//MPI_Barrier(MPI_COMM_WORLD);
	double total = 0;

	if (rank == 0)
	{
		total = S;
		for (int i = 1; i < size; i++)
		{
			MPI_Recv(&S, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			total += S;
		}
		total = sqrt(total);
		printf("E = %f\n", total);
		//printf("E = %f\n", total);
	}
	else MPI_Send(&S, 1, MPI_DOUBLE, 0, 77, MPI_COMM_WORLD);
	//cout << S << endl << endl;
	//MPI_Barrier(MPI_COMM_WORLD);
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