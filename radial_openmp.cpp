#include <iostream>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <new>
#include <string>
#include <algorithm>
#include <sstream>
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
void split_period_double(std::string& input,double& x1,double& x2){
	std::istringstream iss(input);
	std::list<std::string> temp;
	std::string u;
	while(getline(iss,u,' ')){
		temp.push_back(u);
	}
	x1=std::stof(*(temp.begin()));
	x2=std::stof(*(next(temp.begin(),1)));
}
void split_position_double(std::string& input,double& x,double& y,double& z){
	std::istringstream iss(input);
	std::list<std::string> temp;
	std::string u;
	while(getline(iss,u,' ')){
		temp.push_back(u);
	}
	x=std::stof(*(temp.begin()));
	y=std::stof(*(next(temp.begin(),1)));
	z=std::stof(*(next(temp.begin(),2)));
}
void printtraject(int cell,int tick,int threadID_print_signal,std::list<std::string> &input,double variance[][3]){
	  double A[cell*cell*cell][3];
    double B[cell*cell*cell][3];
    double oxygen[3*cell*cell*cell][3];
    std::string la_pattern="ITEM: BOX BOUNDS pp pp pp";
    std::string coord_pattern="ITEM: ATOMS x y z ";
    double* x=new double[3];
    double* period=new double[3];
    double* dist=NULL;
    double sum[3]={0.0,0.0,0.0};
    int* temp_neighbor;
		std::list<double> px;
		std::list<double> py;
		std::list<double> pz;
		int count=0;
    /******************end define those things *************/
    for(std::list<std::string>::iterator line=input.begin();line!=input.end();){
			if(threadID_print_signal==0){
					std::cout<<count++<<std::endl;
			}
			for(size_t i=0;i<3;i++){
							split_period_double(*line,x[0],x[1]);
							line++;
              period[i]=x[1]-x[0];
            }
        for(size_t i=0;i<cell*cell*cell;i++){
					split_position_double(*line,x[0],x[1],x[2]);
					line++;
					for(size_t j=0;j<3;j++){
		    			A[i][j]=x[j];
					}
				}
				for(size_t i=0;i<cell*cell*cell;i++){
					split_position_double(*line,x[0],x[1],x[2]);
					line++;
					for(size_t j=0;j<3;j++){
						  B[i][j]=x[j];
					}
				}
				for(size_t i=0;i<3*cell*cell*cell;i++){
					split_position_double(*line,x[0],x[1],x[2]);
					line++;
					for(size_t j=0;j<3;j++){
						oxygen[i][j]=x[j];
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
    std::fstream dump;
    dump.open("dump.xyz",std::fstream::in);
		int count=0;
		for(std::string line;getline(dump,line);){
			count++;
		}
		dump.close();
    dump.open("dump.xyz",std::fstream::in);
		std::list<std::string> fileall;
		count=0;
    std::string la_pattern="ITEM: BOX BOUNDS pp pp pp";
    std::string coord_pattern="ITEM: ATOMS x y z ";
		for(std::string line;getline(dump,line);){
			if(line==la_pattern){
				for(size_t i=0;i<3;i++){
					getline(dump,line);
					fileall.push_back(line);
				}
			}
			if(line==coord_pattern){
				for(size_t i=0;i<5*cell*cell*cell;i++){
					getline(dump,line);
					fileall.push_back(line);
				}
			}
		}
		std::cout<<"finished reading variables "<<std::endl;
		std::cout<<fileall.size()<<std::endl;
		#pragma omp parallel
		{
			int ID=omp_get_thread_num();
			int tick;
			for(size_t atom_i=0;atom_i<job_per_processor;atom_i++){
				tick=ID*job_per_processor+atom_i;
				printtraject(cell,ID,tick,fileall,variance);
			}
		}
		int tick=0;
		for(size_t atom_i=thread_max*job_per_processor;atom_i<cell*cell*cell;atom_i++){
				tick=atom_i;
				printtraject(cell,0,tick,fileall,variance);		
		}
		std::fstream fs;
		fs.open("final.txt",std::fstream::in);
		for(size_t i=0;i<cell*cell*cell;i++){
			fs<<variance[i][0]<<" "<<variance[i][1]<<" "<<variance[i][2]<<std::endl;
		}
}
