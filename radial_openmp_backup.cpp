#include <iostream>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <new>
#include <string>
#include <omp.h>
#include <list>
int changeback(int nx,int ny,int nz,int p){
     int re;
     re=(nx+p)%p+p*((ny+p)%p)+p*p*((nz+p)%p);
     return re;
}
//compute the distance between two points
double* distance(double* a,double* b,double* p){
	double* dist=new double[3];
	double temp;
	for(size_t i=0;i<3;i++){
		temp=b[i]-a[i];
		temp=(temp/p[i]-round(temp/p[i]))*p[i];
		dist[i]=temp;
	}
	return dist;
}
double average(std::list<double>& input){
		double sum=0; 
			for(std::list<double>::iterator a=input.begin();a!=input.end();a++){
						sum=sum+*a;
							}
				return sum/input.size();
};
double var(std::list<double>& input){
		double sum=0.0;
		double ave=average(input);
		for(std::list<double>::iterator a=input.begin();a!=input.end();a++){
			sum=sum+(*a-ave)*(*a-ave);
		}
		return sum/input.size();
}
//compute the index change from 1D to 3D
int* changeindex(int index,int cell){
	int* re=new int[3];
	re[2]=floor(index/(cell*cell));
	index=index-re[2]*cell*cell;
	re[1]=floor(index/cell);
	re[0]=index-cell*re[1];
	return re;
}
//compute the neighbor oxygen index for A site atoms.
int* neighbor_o_forA(int index,int cell){
	int* index_3D=changeindex(index,cell);
	int* nei=new int[12];
	nei[0]=changeback(index_3D[0],index_3D[1],index_3D[2],cell);
	nei[1]=changeback(index_3D[0]-1,index_3D[1],index_3D[2],cell);
	nei[2]=changeback(index_3D[0],index_3D[1]-1,index_3D[2],cell);
	nei[3]=changeback(index_3D[0]-1,index_3D[1]-1,index_3D[2],cell);
	nei[4]=changeback(index_3D[0],index_3D[1],index_3D[2],cell)+cell*cell*cell;
	nei[5]=changeback(index_3D[0]-1,index_3D[1],index_3D[2],cell)+cell*cell*cell;
	nei[6]=changeback(index_3D[0]-1,index_3D[1],index_3D[2]-1,cell)+cell*cell*cell;
	nei[7]=changeback(index_3D[0],index_3D[1],index_3D[2]-1,cell)+cell*cell*cell;
	nei[8]=changeback(index_3D[0],index_3D[1],index_3D[2],cell)+2*cell*cell*cell;
	nei[9]=changeback(index_3D[0],index_3D[1]-1,index_3D[2],cell)+2*cell*cell*cell;
	nei[10]=changeback(index_3D[0],index_3D[1],index_3D[2]-1,cell)+2*cell*cell*cell;
	nei[11]=changeback(index_3D[0],index_3D[1]-1,index_3D[2]-1,cell)+2*cell*cell*cell;
    return nei;
}
void sum_together(double* sum,double* add,int len){
    for(size_t i=0;i<len;i++){
        sum[i]=sum[i]+add[i];
    }
}
void printtraject(int cell,int tick,int threadIDspeed,std::string input_filename,double variance[][3]){
    double A[cell*cell*cell][3];
    double B[cell*cell*cell][3];
    double oxygen[3*cell*cell*cell][3];
    std::string la_pattern="ITEM: BOX BOUNDS pp pp pp";
    std::string coord_pattern="ITEM: ATOMS x y z ";
    double x1,x2;
    double* period=new double[3];
    double* dist=NULL;
    double sum[3]={0.0,0.0,0.0};
    int* temp_neighbor;
		std::list<double> px;
		std::list<double> py;
		std::list<double> pz;
    /******************end define those things *************/
    std::fstream dump;
		int count=0;
    dump.open(input_filename,std::fstream::in);
    for(std::string line;getline(dump,line);){
      if(la_pattern==line){
				if(threadIDspeed==10){
					std::cout<<count++<<std::endl;
				}
            for(size_t i=0;i<3;i++){
                dump>>x1;
                dump>>x2;
                period[i]=x2-x1;
            }
        }
      if(coord_pattern==line){
        for(size_t i=0;i<cell*cell*cell;i++){
					for(size_t j=0;j<3;j++){
						dump>>A[i][j];
					}
				}
				for(size_t i=0;i<cell*cell*cell;i++){
					for(size_t j=0;j<3;j++){
						dump>>B[i][j];
					}
				}
				for(size_t i=0;i<3*cell*cell*cell;i++){
					for(size_t j=0;j<3;j++){
					    dump>>oxygen[i][j];
					}
      	}
      	temp_neighbor=neighbor_o_forA(tick,cell);
      	for(size_t j=0;j<12;j++){
        	    dist=distance(A[tick],oxygen[temp_neighbor[j]],period);
							sum_together(sum,dist,3);
     		 }
				px.push_back(sum[0]/12);
				py.push_back(sum[1]/12);
				pz.push_back(sum[2]/12);
					for(size_t i=0;i<3;i++){
						sum[i]=0.0;
					}
     		 }
   		 }
		dump.close();
		variance[tick][0]=var(px);
		variance[tick][1]=var(py);
		variance[tick][2]=var(pz);
}
int main(int argc,char* argv[]){
    int cell=20;
		int thread_max=omp_get_max_threads();
		int job_per_processor=floor(cell*cell*cell/thread_max);//how many jobs need to be don for each processor.
		int residue=cell*cell*cell-job_per_processor;//residue job for the main processor.
    double variance[cell*cell*cell][3];
		#pragma omp parallel
		{
	  	std::string input_file=argv[1];
			int ID=omp_get_thread_num();
			int tick;
			for(size_t atom_i=0;atom_i<job_per_processor;atom_i++){
				tick=ID*job_per_processor+atom_i;
				printtraject(cell,tick,ID,argv[1],variance);
			}
		}
		int tick=0;
		for(size_t atom_i=thread_max*job_per_processor;atom_i<cell*cell*cell;atom_i++){
				tick=atom_i;
				printtraject(cell,tick,0,argv[1],variance);		
		}
		std::fstream fs;
		fs.open("final.txt",std::fstream::out);
		for(size_t i=0;i<cell*cell*cell;i++){
			fs<<variance[i][0]<<" "<<variance[i][1]<<" "<<variance[i][2]<<std::endl;
		}
}
