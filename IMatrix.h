/**
  * This virtual class defines the functions that will be supported 
  * by the abstract 2D matric data structure.
  *
  * DO NOT CHANGE THIS CLASS!
  */
class IMatrix
{
	public:
    
	/**
	  * This initializes an array that will contain N rows and M columns. 
	  * The values in the initial matrix are assumed to be 0. 
	  * This method must be called exactly once when creating the matrix, and 
	  * subsequent calls to this method will should return a "logic_error" exception. 
	  */
    virtual void init (int , int )=0; 

	/**
	  * If append is called for the first time, it sets the values of row 0 to the 
	  * values passed in row_vals. Else, it sets the values of the next row of the 
	  * matrix since the last call to append to these values. row_vals must have 
	  * exactly M items in it, and you cannot append more then N rows to the matrix. 
	  * If either of these are violated, you should throw an "out_of_range" exception.
	  */
    virtual void append (float*)=0; 
	
	/**
	  * This returns the value present in row i and column j of the matrix, 
	  * indexed from 0. You should throw an "out_of_range" exception if the row or 
	  * column values are outside the valid range of rows and columns.
	  */
    virtual float get(int , int )=0; 

	/**
	  * This returns the number of rows in the matrix based on init()
	  */
	virtual int row_count()=0;
	
	/**
	  * This returns the number of columns in the matrix based on init()
	  */
	virtual int col_count()=0;
	
	/**
	  * This returns the number of non zero values in the matrix
	  */
	virtual int nnz_count()=0;
	
};

