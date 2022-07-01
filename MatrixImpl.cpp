#include"IMatrix.h"
#include<iostream>
#include<stdexcept>
#include<vector>
#include<cstdlib>

using namespace std;


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * 2D Matrix Implementation usign Arrays
 */
class ArrayMatrix: public IMatrix
{  
	/*#######################################################
	 # Define your own data structures here 
	 #######################################################*/

	private: 
		// TODO
		float **mat;
		int nrows, ncols, nnz, append_cnt;
   
	public:
		// TODO
   
   
	/**--
     * Default constructor is required.
	 * You may optionally define additional constructors as well.
     */
	ArrayMatrix()
	{ 
		// TODO: Provide implementation for default constructor
		mat = NULL;
		nnz = 0;
		append_cnt = 0;
		nrows = 0;
		ncols = 0;
	}

 
 
	/*#######################################################
 	  # Override functions defined in IMatrix interface
	  #######################################################*/

	/** Override the init function defined in the IMatrix interface */
	void init (int N, int M)
 	{ 
		// TODO: Provide implementation for init using array

		if(mat != NULL) throw logic_error("Matrix can be initialized only once.");

		nrows = N;
		ncols = M;

 		mat = (float**) malloc(nrows*sizeof(float*));
		for(int i=0; i<nrows; i++)
		{
			//calloc allocates memory and sets bits to 0
			mat[i] = (float*) calloc(ncols, sizeof(float));
		}

	}

	
	/** Override the append function defined in the IMatrix interface */
   void append (float* row_vals)
	{ 
	    // TODO: Provide implementation for append using array
		for(int i=0; i<ncols; i++)
		{
			if(row_vals[i]) nnz++;
			mat[append_cnt][i] = row_vals[i];
		}
		append_cnt++;
	}

	
	/** Override the get function defined in the IMatrix interface */
	float get(int i, int j)
	{ 
	    // TODO: Provide implementation for get using array
	    if(i<0 || i>=nrows || j<0 || j>=ncols) throw out_of_range("Array index out of range");
	    return mat[i][j];
	}

	
	/**
	  * This returns the number of rows in the matrix based on init()
	  */
	int row_count()
	{
		// TODO: Provide implementation using array
		return nrows;
	}
	
	/**
	  * This returns the number of columns in the matrix based on init()
	  */
	int col_count() 
	{
		// TODO: Provide implementation using array
		return ncols;
	}

	
	/**
	  * This returns the number of non zero values in the matrix
	  */
	int nnz_count()
	{
		// TODO: Provide implementation using array
		return nnz;
	}

	/*#######################################################
	 # Optionally define any other private utility functions here
	 #######################################################*/

	// TODO
};

/////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * 2D Matrix Implementation usign CSR
 */
class CSRMatrix: public IMatrix
{  	

	/*#######################################################
	 # Define your own data structures here 
	 #######################################################*/

	private: 
		// TODO
		vector<float> a;
		vector<int> col_id;
		int *row_nnz_count;
		int nrows, ncols, append_cnt;
   
	public:
		// TODO
   
   
	/**
     * Default constructor is required.
	 * You may optionally define additional constructors as well.
     */
	CSRMatrix()
	{ 
		// TODO: Provide implementation for default constructor
		row_nnz_count = NULL;
 		append_cnt = 0;
 		nrows = 0;
 		ncols = 0;
	}
 	
 
	/*#######################################################
 	  # Override functions defined in IMatrix interface
	  #######################################################*/

	/** Override the init function defined in the IMatrix interface */
	void init (int N, int M)
 	{ 
		// TODO: Provide implementation for init using csr
		if(row_nnz_count != NULL) throw logic_error("Matrix can be initialized only once.");

 		row_nnz_count = (int *)malloc((N+1)*sizeof(int));
 		nrows = N;
 		ncols = M;
 		row_nnz_count[0] = 0;
	}
	

	/** Override the append function defined in the IMatrix interface */
   void append (float* row_vals)
	{ 
	    // TODO: Provide implementation for append using csr
	    int nnz=0;
	    for(int i=0; i<ncols; i++)
	    {
	    	// add to a only if element is non zero
	    	if(row_vals[i])
	    	{
	    		a.push_back(row_vals[i]);
	    		col_id.push_back(i);
	    		nnz++;
	    	}    		
	    }
	    append_cnt++;
	    row_nnz_count[append_cnt] = row_nnz_count[append_cnt-1]+nnz;
	}

	
	/** Override the get function defined in the IMatrix interface */
	float get(int i, int j)
	{ 
	    // TODO: Provide implementation for get using csr
	    if(i<0 || i>=nrows || j<0 || j>=ncols) throw out_of_range("Array index out of range");

	    int id = get_id(i,j);

	    return (id >= 0)? a[id] : 0.0;
	}

	/**
	  * This returns the number of rows in the matrix based on init()
	  */
	int row_count()
	{
		// TODO: Provide implementation using csr
		return nrows;
	}
	
	/**
	  * This returns the number of columns in the matrix based on init()
	  */
	int col_count() 
	{
		// TODO: Provide implementation using csr
		return ncols;
	}

	
	/**
	  * This returns the number of non zero values in the matrix
	  */
	int nnz_count()
	{
		// TODO: Provide implementation using csr
		return row_nnz_count[nrows];
	}

	
	/*#######################################################
	 # Optionally define any other private utility functions here
	 #######################################################*/
	
	// TODO
	private:
	
	// this function returns index in a corresponding to given pair (i,j)
	int get_id(int i, int j)
	{
		// number of non zero elements before row i
		int offset = row_nnz_count[i];

		// number of non zero elements in current row i
		int nnz_cnt = row_nnz_count[i+1] - offset;

		// return the index in array a if column matches else return a -1;
		while(nnz_cnt--)
		{
			if(col_id[offset++] == j) return offset-1;
		}

		return -1;
	}
	
};
