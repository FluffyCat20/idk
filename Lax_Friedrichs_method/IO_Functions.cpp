/*//#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdarg.h> //для va_list
//#include "headers.h"
#include <float.h>	//для DBL_MAX*/

/*#ifdef __linux //for linux compiler
#include <sys/stat.h>
#include <sys/types.h>
#define __int32 int32_t
#define __int64 int64_t
#elif _WIN32 //for windows compiler
#include <windows.h>
#endif*/

/*#define TECPLOT_STRING_BUFFER_LENGTH 64 //длина буфера для записи строк в бинарном формате
#define TECPLOT_FIELDS_NUMBER 2 + INT_EQUATIONS_COUNT + 1 //число величин для записи в tecplot: 2 координаты, все величины + число Маха*/

#include "IO_Functions.h"

//binery writing functions for .plt files
int writeBinary4BytesInt(__int32 value, FILE* file)
{
	__int32 dword;
	dword = value;
	fwrite(&dword, 4, 1, file);
	return 0;
}
int writeBinary4BytesFloat(float value, FILE* file)
{
	__int32 dword;
	*((float*)&dword) = value;
	fwrite(&dword, 4, 1, file);
	return 0;
}
int writeBinary8BytesDouble(double value, FILE* file)
{
	__int64 qword;
	*((double*)&qword) = value;
	fwrite(&qword, 8, 1, file);
	return 0;
}
int writeBinaryString(char* str, FILE* file)
{
	__int32 dword;
	int p;
	p = 0;
	do
	{
		dword = str[p++];
		fwrite(&dword, 4, 1, file);
	} while (dword);
	return 0;
}

void WriteAsTecplotFile(char* filename)
{
	FILE* fp;
	fp = fopen(filename, "wb");
	WriteTecplotFileHeader(fp);
	WriteTecplotZoneHeader(fp);
	WriteTecplotFinFileHeader(fp);
	WriteTecplotZoneDataFloat(fp);
	fclose(fp);
}

void WriteTecplotFileHeader(FILE *fp)
{
	char str[TECPLOT_STRING_BUFFER_LENGTH] = "";
	fwrite("#!TDV111", 8, 1, fp); //Код версии данных Tecplot, нужно именно "1 символ = 1 байт"
	writeBinary4BytesInt(1, fp); //byte order
	writeBinary4BytesInt(0, fp); //filetype: 0 = full, 1 = grid, 2 = solution
	strcpy(str, "Frame title"); //frame title (visible if "frame header" is enabled)
	writeBinaryString(str, fp);
	writeBinary4BytesInt(TECPLOT_FIELDS_NUMBER, fp); //Число массивов, включая координаты
	strcpy(str, "x"); writeBinaryString(str, fp); //Названия параметров
	strcpy(str, "y"); writeBinaryString(str, fp);
	strcpy(str, "<greek>r</greek>"); writeBinaryString(str, fp);
	strcpy(str, "U"); writeBinaryString(str, fp);
	strcpy(str, "V"); writeBinaryString(str, fp);
	strcpy(str, "P"); writeBinaryString(str, fp);
	strcpy(str, "M"); writeBinaryString(str, fp);
}

void WriteTecplotZoneHeader(FILE *fp)
{
	__int32 tmp;
	char str[TECPLOT_STRING_BUFFER_LENGTH] = "";
	writeBinary4BytesFloat(299.0f, fp); //zone marker
	sprintf(str, "%.10f (%d)", double_time, int_my_process_id); //название зоны
	writeBinaryString(str, fp);
	writeBinary4BytesInt(-1, fp); //parent zone
	writeBinary4BytesInt(-2, fp); //time strand (-2 = auto, -1 = none, 0+ = explicit value)
	writeBinary8BytesDouble(double_time, fp); //время в формате double
	writeBinary4BytesInt(-1, fp); //not used, set to -1
	writeBinary4BytesInt(0, fp); //zone type, 0 = ordered
	writeBinary4BytesInt(0, fp); //0 - по блокам, 1 - по точкам (по точкам в бинарном похоже не поддерживатеся)
	writeBinary4BytesInt(0, fp); //var location, 0 = all data is at the nodes
	writeBinary4BytesInt(0, fp); //0 = no face neighbor data supplied
	writeBinary4BytesInt(0, fp); //0 = no user defined face connections
	writeBinary4BytesInt(1, fp); //число узлов вдоль оси Z (Imax), fastest index
	writeBinary4BytesInt(int_block_Lmax + 1 - 2*int_weno_r, fp); //узлы вдоль оси Y (Jmax), slower index
	writeBinary4BytesInt(int_block_Kmax + 1 - 2*int_weno_r, fp); //узлы вдоль оси X (Kmax), slowest index
	writeBinary4BytesInt(0, fp); //no auxiliary names/values
}

void WriteTecplotFinFileHeader(FILE *fp)
{
	writeBinary4BytesFloat(357.0f, fp); //end of header marker
}

void WriteTecplotZoneDataFloat(FILE *fp)
{
	int k, l, p;
	__int32 dword;
	__int64 qword;
	char str[TECPLOT_STRING_BUFFER_LENGTH] = "";
	double minval, maxval, double_tmp;
	writeBinary4BytesFloat(299.0f, fp); //zone marker
	//коды типов данных для каждой величины
	for (p = 0; p < TECPLOT_FIELDS_NUMBER; p++)
		writeBinary4BytesInt(1, fp); //1 = float
	writeBinary4BytesInt(0, fp); //no passive variables
	writeBinary4BytesInt(0, fp); //no variable sharing
	writeBinary4BytesInt(-1, fp); //no connectivity list zone
	//минимумы и максимумы
	writeBinary8BytesDouble(double_block_x_min, fp);
	writeBinary8BytesDouble(double_block_x_max, fp);
	writeBinary8BytesDouble(double_block_y_min, fp);
	writeBinary8BytesDouble(double_block_y_max, fp);
	for (p = 2; p < TECPLOT_FIELDS_NUMBER; p++)
	{
		minval = DBL_MAX;
		maxval = -DBL_MAX;
		for (k = int_weno_r; k <= int_block_Kmax - int_weno_r; k++)
		for (l = int_weno_r; l <= int_block_Lmax - int_weno_r; l++)
		{
			switch (p)
			{
			case 2:
				minval = R[k][l] < minval ? R[k][l] : minval;
				maxval = R[k][l] > maxval ? R[k][l] : maxval;
				break;
			case 3:
				minval = U[k][l] < minval ? U[k][l] : minval;
				maxval = U[k][l] > maxval ? U[k][l] : maxval;
				break;
			case 4:
				minval = V[k][l] < minval ? V[k][l] : minval;
				maxval = V[k][l] > maxval ? V[k][l] : maxval;
				break;
			case 5:
				minval = P[k][l] < minval ? P[k][l] : minval;
				maxval = P[k][l] > maxval ? P[k][l] : maxval;
				break;
			case 6:
				double_tmp = sqrt(U[k][l] * U[k][l] + V[k][l] * V[k][l]) / sqrt(double_gamma*P[k][l] / R[k][l]);
				minval = double_tmp < minval ? double_tmp : minval;
				maxval = double_tmp > maxval ? double_tmp : maxval;
				break;
			}
		}
		writeBinary8BytesDouble(minval, fp);
		writeBinary8BytesDouble(maxval, fp);

	}
	//массивы значений
	for (p = 0; p < TECPLOT_FIELDS_NUMBER; p++)
	{
		for (k = int_weno_r; k <= int_block_Kmax - int_weno_r;)
		{
			for (l = int_weno_r; l <= int_block_Lmax - int_weno_r;)
			{
				switch (p)
				{
				case 0:
					writeBinary4BytesFloat(double_block_x_min + (k - int_weno_r) * double_grid_step, fp);
					break;
				case 1:
					writeBinary4BytesFloat(double_block_y_min + (l - int_weno_r) * double_grid_step, fp);
					break;
				case 2:
					writeBinary4BytesFloat(R[k][l], fp);
					break;
				case 3:
					writeBinary4BytesFloat(U[k][l], fp);
					break;
				case 4:
					writeBinary4BytesFloat(V[k][l], fp);
					break;
				case 5:
					writeBinary4BytesFloat(P[k][l], fp);
					break;
				case 6:
					writeBinary4BytesFloat(sqrt(U[k][l] * U[k][l] + V[k][l] * V[k][l]) / sqrt(double_gamma * P[k][l] / R[k][l]), fp);
					break;
				}
				if (int_export_grid_cell_skip > 1 && l < int_block_Lmax - int_weno_r)
				if (l + int_export_grid_cell_skip > int_block_Lmax - int_weno_r)
					l = int_block_Lmax - int_weno_r;
				else
					l += int_export_grid_cell_skip;
				else
					l++;
			}
			if (int_export_grid_cell_skip > 1 && k < int_block_Kmax - int_weno_r)
			if (k + int_export_grid_cell_skip > int_block_Kmax - int_weno_r)
				k = int_block_Kmax - int_weno_r;
			else
				k += int_export_grid_cell_skip;
			else
				k++;
		}
	}
}
