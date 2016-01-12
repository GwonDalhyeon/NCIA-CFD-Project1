#pragma once

#include <iostream>
#include <cmath>
#include <fstream>

#include "Array2D.h"
#include "VectorND.h"

using namespace std;
template <class TT>
class CSR
{
public:
	//Array2D<TT> values;
	//Array2D<int> columns;
	//Array2D<int> indPrt;

	VectorND<TT> values;
	VectorND<int> columns;
	VectorND<int> indPrt;
	
	int colNum;
	int rowNum;
	int valueNum;
	
	CSR();
	~CSR();

	CSR(const Array2D<TT>& ipArray);

	inline void operator = (const CSR<TT>& ipCSR)
	{
		values = ipCSR.values;
		columns = ipCSR.columns;
		indPrt = ipCSR.indPrt;

		colNum = ipCSR.colNum;
		rowNum = ipCSR.rowNum;
		valueNum = ipCSR.valueNum;
	}

private:

};
template <class TT>
CSR<TT>::CSR()
{
}

template <class TT>
CSR<TT>::~CSR()
{
	cout << "delete CSR" << endl;
}

template<class TT>
inline CSR<TT>::CSR(const Array2D<TT>& ipArray)
{
	rowNum = ipArray.iRes;
	colNum = ipArray.jRes;
	indPrt = VectorND<int>(rowNum);

	TT* tempVal = new double[int(floor(sqrt(double(rowNum*colNum)))) * 10];
	int* tempCol = new int[int(floor(sqrt(double(rowNum*colNum)))) * 10];
	

	for (int i = 0; i < rowNum + 1; i++)
	{
		indPrt.values[i] = -1;

	}
	int tempIndex = 0;


	for (int i = 0; i < rowNum; i++)
	{
		for (int j = 0; j < colNum; j++)
		{
			if (ipArray(i,j) != 0)
			{
				tempVal[tempIndex] = ipArray(i, j);
				tempCol[tempIndex] = j;
				if (indPrt.values[i]<0)
				{
					indPrt.values[i] = tempIndex;
				}
				tempIndex = tempIndex + 1;
			}
		}
	}
	valueNum = tempIndex;
	//indPrt.values[rowNum] = tempIndex;

	values = VectorND<TT>(valueNum);
	columns = VectorND<int>(valueNum);

	for (int i = 0; i < valueNum; i++)
	{
		values.values[i] = tempVal[i];
		columns.values[i] = tempCol[i];
	}

	delete[] tempVal, tempCol;
}

template<class TT>
inline std::ostream& operator<<(std::ostream& output, const CSR<TT>& ipCSR)
{
	output << "CSR" << endl;
	output << "- colNum = " << ipCSR.colNum << endl;
	output << "- rowNum = " << ipCSR.rowNum << endl;
	output << "- valueNum = " << ipCSR.valueNum << endl;

	return output;
}
