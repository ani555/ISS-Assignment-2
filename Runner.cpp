#include <iostream>
#include <string.h>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <stdexcept>
#include <time.h>
#include <list>
#include <queue>
#include <utility>
#include <fstream>
#include <sstream>
#include <algorithm>
#include"MatrixImpl.cpp"

using namespace std;

void add(char* mat_type, char* input_file1, char* input_file2, char* output_file);
void add(IMatrix* mat1, IMatrix* mat2, IMatrix* mat3);
void multiply(IMatrix* mat1, IMatrix* mat2, IMatrix* mat3);
void multiply(char* mat_type, char* input_file1, char* input_file2, char* output_file);
/**
  * Runner code that serves as harness for invoking various functions required 
  * to be implemented for this assignment.
  * You may modify this code, but need to retain the input parameter signature and 
  * output requirements specified in the assignment.
  */  
  
/////////////////////////////////////////////////////////////////////////////
// helpers 
int parseLine(char *line)
{
	int i = strlen(line);
	const char*p = line;
	while(*p < '0' || *p > '9') p++;
	line[i-3] = '\0';
	i = atoi(p);
	return i;
}
//getValue() is a helper function to calculate the RAM usage by the current process: source stackoverflow
int getValue()
{
	FILE *file = fopen("/proc/self/status","r");
	if(file==NULL) perror("Couldnot open file");
	int result = -1;
	char line[128];

	while(fgets(line, 128, file) != NULL)
	{
		if(strncmp(line,"VmRSS:",6) == 0)
		{
			result = parseLine(line);
			break;
		}
	}
	fclose(file);
	return result;
}

void append_to_file(char *file_name, int sz, double val)
{
	FILE *file = fopen(file_name,"a");
	fprintf(file, "%d,%.4lf\n", sz, val);
	fclose(file);
}

void append_to_file(char *file_name, int sz, int val)
{
	FILE *file = fopen(file_name,"a");
	fprintf(file, "%d,%d\n", sz, val);
	fclose(file);
}

void print(IMatrix *mat)
{
	int N = mat->row_count();
	int M = mat->col_count();
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<M; j++)
		{
			printf("%f\t",mat->get(i,j));
		}
		printf("\n");
	}
}

// main code starts here

struct timespec start, finish;


IMatrix* load_arr(char* input_file) 
{
	// TODO: Code for loading matrix from input_file into array impl and returning object
	
	int N , M;
	ifstream in_file(input_file);
	string line;
	int rcnt = 0;
	ArrayMatrix *mat = new ArrayMatrix();
	float *curr_row;
	int dims[2];
	while(getline(in_file, line))
	{
		istringstream iss(line);
		string token;
		int ntokens = 0;
		// The first line should contain only positive integers separated by tab
		// And the number of tokes should be exactly 2
		// All the remaining lines should contain exactly M tokens and tokens must be valid floating pts
		// There must be total N+1 rows
		if(rcnt == 0)
		{
			while(getline(iss, token, '\t'))
			{
				ntokens++;
				//check for valid token
				if(token.find_first_not_of("0123456789") != string::npos) throw invalid_argument("Not a valid matrix dimension");
				int num = atoi(token.c_str());
				if(ntokens > 2) throw invalid_argument("Expected two dimensions");
				dims[ntokens-1] = num;	
			}
			if(ntokens < 2) throw invalid_argument("Expected two dimensions");
			N = dims[0];
			M = dims[1];
			mat->init(N,M);
			curr_row = new float[M];
		}
		else
		{
			while(getline(iss, token, '\t'))
			{
				ntokens++;
				//check for valid token
				if(token.find_first_not_of("-.0123456789") != string::npos) throw invalid_argument("Not a valid float type");
				float num = atof(token.c_str());
				if(ntokens > M) throw invalid_argument("Expected M floats in a row");
				curr_row[ntokens-1] = num;
			}
			if(ntokens < M) throw invalid_argument("Expected M floats in a row");
			mat->append(curr_row);
		}
		rcnt++;
	}

	if(rcnt != N+1) throw invalid_argument("Expected N+1 total rows");
	in_file.close();

	return (IMatrix *) mat;
}

IMatrix* load_csr(char* input_file) 
{

	// TODO: Code for loading matrix from input_file into CSR impl and returning object

	int N , M;
	ifstream in_file(input_file);
	string line;
	int rcnt = 0;
	CSRMatrix *mat = new CSRMatrix();
	float *curr_row;
	int dims[2];
	while(getline(in_file, line))
	{
		istringstream iss(line);
		string token;
		int ntokens = 0;
		// The first line should contain only positive integers separated by tab
		// And the number of tokes should be exactly 2
		// All the remaining lines should contain exactly M tokens and tokens must be valid floating pts
		// There must be total N+1 rows
		if(rcnt == 0)
		{
			while(getline(iss, token, '\t'))
			{
				ntokens++;
				//check for valid token
				if(token.find_first_not_of("0123456789") != string::npos) throw invalid_argument("Not a valid matrix dimension");
				int num = atoi(token.c_str());
				if(ntokens > 2) throw invalid_argument("Expected two dimensions");
				dims[ntokens-1] = num;	
			}
			if(ntokens < 2) throw invalid_argument("Expected two dimensions");
			N = dims[0];
			M = dims[1];
			mat->init(N,M);
			curr_row = new float[M];
		}
		else
		{
			while(getline(iss, token, '\t'))
			{
				ntokens++;
				//check for valid token
				if(token.find_first_not_of("-.0123456789") != string::npos) throw invalid_argument("Not a valid float type");
				float num = atof(token.c_str());
				if(ntokens>M) throw invalid_argument("Expected M floats in a row");
				curr_row[ntokens-1] = num;
			}
			if(ntokens < M) throw invalid_argument("Expected M floats in a row");
			mat->append(curr_row);
		}
		rcnt++;
	}
	in_file.close();
	if(rcnt != N+1) throw invalid_argument("Expected N+1 total rows");
	return (IMatrix *) mat;
}

IMatrix* init_arr(int rows, int cols) 
{
	// TODO: Code for initializing an empty matrix using array impl with rows and cols as 
	// dimensions, and returning the object
	ArrayMatrix *arr_mat = new ArrayMatrix();
	arr_mat->init(rows, cols);
	return (IMatrix *)arr_mat;
}


IMatrix* init_csr(int rows, int cols) 
{
	// TODO: Code for initializing an empty matrix using CSR impl with rows and cols as 
	// dimensions, and returning the object
	CSRMatrix *csr_mat = new CSRMatrix();
	csr_mat->init(rows, cols);
	return (IMatrix *)csr_mat;
}


void print_mat(IMatrix* mat, char* output_file) 
{
	// TODO: print matrix as TSV to otput_file
	FILE *file = fopen(output_file, "w");
	int N = mat->row_count();
	int M = mat->col_count();
	
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<M; j++)
		{
			if(j==M-1)
				fprintf(file, "%.3f\n",mat->get(i,j));
			else
				fprintf(file, "%.3f\t",mat->get(i,j));
		}
	}
	fclose(file);

}


/////////////////////////////////////////////////////////////////////////////

void load(char* mat_type, char* input_file, char* output_file)
{
	// TODO: any other validation?

	if(input_file==NULL|| output_file==NULL) throw invalid_argument("File passed is NULL.");
	if(strcmp(".tsv", input_file+(strlen(input_file)-4)) ) throw invalid_argument("Input file passed is invalid.");
	if(strcmp(".tsv", output_file+(strlen(output_file)-4)) ) throw invalid_argument("Output file passed is invalid.");
	
	IMatrix* mat1;

	if (strcmp("array", mat_type)==0) // TODO: time this region and print "load,array,output_file,time_millisec"
    {
    	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    	mat1 = load_arr(input_file);
    	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &finish);
    	double run_time = ((finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/1e9)*1e3;
    	printf("load,array,%s,%.4lf\n",output_file,run_time);
    	
		if(mat1!=NULL)
		{
			// Print the memory
			char r[30] = "op/load_array_mem.csv";
			append_to_file(r, mat1->row_count(), getValue());
		
			// Print the time
			char s[30] = "op/load_array_time.csv";
    		append_to_file(s, mat1->row_count(), run_time);
    	}

    }  
	else if (strcmp("csr", mat_type)==0) // TODO: time this region and print "load,csr,output_file,time_millisec"
	{  
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
		mat1 = load_csr(input_file);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &finish);
		double run_time = ((finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/1e9)*1e3;
		printf("load,csr,%s,%0.4lf\n",output_file,run_time);	

		if(mat1!=NULL)
		{
			// Print the memory
			char r[30] = "op/load_csr_mem.csv";
			append_to_file(r, mat1->row_count(), getValue());
		
			// Print the time
			char s[30] = "op/load_csr_time.csv";
    		append_to_file(s, mat1->row_count(), run_time);
    	}

	}
	else
		cout<<"[load] invalid matrix type input seen: "<<mat_type<<endl;

	// store the loaded matrix mat1 in output file given by output_file
	print_mat(mat1, output_file);


	return;
}

/////////////////////////////////////////////////////////////////////////////

void add(char* mat_type, char* input_file1, char* input_file2, char* output_file)
{ 
	// TODO: any other validation?
	if(input_file1==NULL || input_file2==NULL || output_file==NULL)
		throw invalid_argument("File passed is NULL");


	IMatrix *mat1, *mat2, *mat3;
	
	if (strcmp("array", mat_type)==0) 
    {
		mat1 = load_arr(input_file1);
		mat2 = load_arr(input_file2);
		
		// TODO: any other validation?
		if(mat1->row_count()!=mat2->row_count() || mat1->col_count()!=mat2->col_count())
			throw logic_error("Dimensions of the input matrices dont match");

		mat3 = init_arr(mat1->row_count(), mat1->col_count());
    }  
	else if (strcmp("csr", mat_type)==0)
	{
		mat1 = load_csr(input_file1);
		mat2 = load_csr(input_file2);

		// TODO: any other validation?
		if(mat1->row_count()!=mat2->row_count() || mat1->col_count()!=mat2->col_count())
			throw logic_error("Dimensions of the input matrices dont match");
		
		mat3 = init_csr(mat1->row_count(), mat1->col_count());
	}
	else {
		cout<<"[add] invalid matrix type input seen: "<<mat_type<<endl;
		return;
	}
	
	// TODO: time this method and print "add,mat_type,output_file,time_millisec"
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);	
	add(mat1, mat2, mat3);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &finish);
	double run_time = ((finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/1e9)*1e3;
	printf("add,%s,%s,%0.4lf\n",mat_type,output_file,run_time);			
	
	// store output matrix mat3 in file given by output_file
	print_mat(mat3, output_file);

	return;
}

void add(IMatrix* mat1, IMatrix* mat2, IMatrix* mat3) 
{	
	// TODO: Code for adding the mat1 and mat2 and storing in a third matrix mat3
	int N = mat1->row_count();
	int M = mat1->col_count();

	float *curr_row = (float *)malloc(M*sizeof(float));
	for(int i=0; i<N; i++)
	{		
		for(int j=0; j<M; j++)
		{
			curr_row[j] = mat1->get(i,j) + mat2->get(i,j);
		}
		mat3->append(curr_row);
	}
	return;
}

/////////////////////////////////////////////////////////////////////////////

void multiply(char* mat_type, char* input_file1, char* input_file2, char* output_file)
{
	// TODO: any other validation?
	if(input_file1==NULL || input_file2==NULL || output_file==NULL)
		throw invalid_argument("File passed is NULL");

	IMatrix *mat1, *mat2, *mat3;

	if (strcmp("array", mat_type)==0) 
    { 
		mat1 = load_arr(input_file1);
		mat2 = load_arr(input_file2);
		
		// TODO: any other validation?
		if(mat1->col_count()!=mat2->row_count())
			throw logic_error("Invalid dimensions for matrix multiplication");
		
		// TODO: init output matrix mat3 with arr impl
		mat3 = init_arr(mat1->row_count(), mat2->col_count());
	}  
	else if (strcmp("csr", mat_type)==0)
	{
		mat1 = load_csr(input_file1);
		mat2 = load_csr(input_file2);

		// TODO: any other validation?
		if(mat1->col_count()!=mat2->row_count())
			throw logic_error("Invalid dimensions for matrix multiplication");

		
		// TODO: init output matrix mat3 with csr impl
		mat3 = init_csr(mat1->row_count(), mat2->col_count());
	}
	else {
		cout<<"[multiply] invalid matrix type input seen: "<<mat_type<<endl;
		return;
	}
	
	// TODO: time this method and print "multiply,csr,output_file,time_millisec"
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);	
	multiply(mat1, mat2, mat3);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &finish);
	double run_time = ((finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/1e9)*1e3;
	printf("multiply,%s,%s,%.4lf\n",mat_type,output_file,run_time);	

		
	// store output matrix mat3 in file given by output_file
	print_mat(mat3, output_file);

	return;
}

void multiply(IMatrix* mat1, IMatrix* mat2, IMatrix* mat3) 
{	
	
	// TODO: Code for adding the mat1 and mat2 and storing in a third matrix mat3
	// If mat1 has dimension NxM and mat2 has MxP mat3 has dimensions NxP
	int N = mat1->row_count();
	int M = mat1->col_count();
	int P = mat2->col_count();

	float *curr_row = (float *) malloc(P*sizeof(float));

	for(int i=0; i<N; i++)
	{
		for(int j=0; j<P; j++)
		{
			float tmp = 0.0;
			for(int k=0; k<M; k++)
			{
				tmp+= mat1->get(i,k)*mat2->get(k,j);
			}
			curr_row[j] = tmp;
		}
		mat3->append(curr_row);
	}

	return;
}
/////////////////////////////////////////////////////////////////////////////

void bfs(char* input_file, char* root_id, char* output_file)
{
 
	// TODO: any validation?
	if(input_file==NULL || output_file==NULL)
		throw invalid_argument("File passed is NULL");

	if(root_id==NULL)
		throw invalid_argument("Root cannot be NULL");

	IMatrix* mat1;
	// TODO: Define a List ADT traverse_list to store output.
	// TODO
	mat1 = load_csr(input_file);

	if(mat1->row_count()!=mat1->col_count())
		throw logic_error("Row and column dimensions of adjacency matrix not equal");

	list<pair<int, list<int> > > depth_vids;
	queue<int> q;
	int N = mat1->row_count();
	bool visited[N];
	int depth[N];
	// Do the initializations
	memset(visited, false, sizeof(visited));
	int cols = mat1->col_count();
	int root = atoi(root_id);
	q.push(root);
	depth[root]=0;
	visited[root]=true;

	// TODO: time this region and print "bfs,csr,output_file,time_millisec"
	// TODO: Code for doing BFS on the matrix starting with vertex present in row "root_id"
	// TODO: Add traversed items into traverse_list
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
	while(!q.empty())
	{	
		
		int v = q.front();
		q.pop();
		//if already that depth is present then add the present vertex to that list
		list<pair<int,list<int> > >::iterator it;
		bool updated = false;
		for(it=depth_vids.begin(); it!=depth_vids.end(); ++it)
		{
			if((*it).first==depth[v])
			{
				(*it).second.push_back(v);
				updated=true;
				break;
			}
		}
		// if depth not present create a new entry in depth vids with list<vids> containing the present vertex
		if(!updated)
		{
			pair<int, list<int> > new_v;
			new_v.first = depth[v];
			new_v.second.push_back(v);
			depth_vids.push_back(new_v);	
		}

		// visit neighbours of v and push them to q
		

		for(int i=0; i<cols; i++)
		{
			if(mat1->get(v, i) == 1.0 && !visited[i])
			{
				visited[i] = true;
				q.push(i);
				depth[i] = depth[v] + 1;
			}
		}
	}
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &finish);
	double run_time = ((finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/1e9)*1e3;
	printf("bfs,csr,%s,%0.4lf\n",output_file,run_time);	

	
	// TODO: store traversal output present in traverse_list into file given by output_file
	list<pair<int, list<int> > >:: iterator it1;
	list<int>:: iterator it2;
	FILE *file = fopen(output_file, "w");
	for(it1 = depth_vids.begin(); it1 != depth_vids.end(); it1++)
	{
		int d = (*it1).first;
		list<int> vids = (*it1).second;
		vids.sort();
		int vertex_cnt = 0;
		fprintf(file,"%d,",d);
		for(it2 = vids.begin(); it2 != vids.end(); it2++)
		{
			vertex_cnt++;
			if(vertex_cnt == vids.size())
				fprintf(file,"%d\n",*it2);
			else
				fprintf(file,"%d,",*it2);
 
		}
	}
	return;

}


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

int main(int n, char *argv[])
{

	if(strcmp("load", argv[1])==0)
	{
		load(argv[2], argv[3], argv[4]);
	}
    else if( strcmp("add", argv[1])==0)
	{
		add(argv[2], argv[3], argv[4], argv[5]);
	}
    else if( strcmp("multiply", argv[1])==0 )
	{
        multiply(argv[2], argv[3], argv[4], argv[5]);
	}
    else if(strcmp("bfs", argv[1])==0)
	{
        bfs(argv[2], argv[3], argv[4]);
	} else 
		cout<<"[main] invalid input parameters. Valid usage is..."<<endl;


	return 0;
    
}

