#pragma once

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdarg.h> //для va_list
#include <float.h>	//для DBL_MAX

#define TECPLOT_STRING_BUFFER_LENGTH 64 //длина буфера для записи строк в бинарном формате
#define TECPLOT_FIELDS_NUMBER 2 + 4 + 1 //число величин для записи в tecplot: 2 координаты, все величины + число Маха

int writeBinary4BytesInt(__int32 value, FILE* file);

int writeBinary4BytesFloat(float value, FILE* file);

int writeBinary8BytesDouble(double value, FILE* file);

int writeBinaryString(char* str, FILE* file);

void WriteAsTecplotFile(char* filename);

void WriteTecplotFileHeader(FILE *fp);

void WriteTecplotZoneHeader(FILE *fp);

void WriteTecplotFinFileHeader(FILE *fp);

void WriteTecplotZoneDataFloat(FILE *fp);
