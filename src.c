#include<stdio.h>
#include<mpi.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

int main( int argc, char *argv[])
{
int fun=atoi (argv[3]);

if(fun==1)
{
int myrank, size,N,rw,c,pr,T;
int p;
double *db,*rb,*lb,*ub;

  MPI_Status status;
double sTime,eTime,time,max=0.0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  
  p=size;
  N=atoi (argv[1]);
  T=atoi (argv[2]);
  
  rw=sqrt(p);
  pr=sqrt(N);
 
 double **mat;
 
 mat=malloc(pr*sizeof(double*));
 mat[0]=malloc(pr*pr*sizeof(double));
 
 for(int i=1;i<pr;i++)
 mat[i]=mat[0]+pr*i;
 
 
 for(int i=0;i<pr;i++)
 {
     for(int j=0;j<pr;j++)
     {
         mat[i][j]=rand()%10;
     }
 }
 

lb= (double *) malloc(pr * sizeof(double));   
rb=(double *) malloc(pr * sizeof(double));
ub=(double *) malloc(pr * sizeof(double));
db=(double *) malloc(pr * sizeof(double));

 memset(lb, 0, pr * sizeof(double));
 memset(rb, 0, pr * sizeof(double));
 memset(ub, 0, pr * sizeof(double));
 memset(db, 0, pr * sizeof(double));

 
MPI_Request reqs[N];

sTime=MPI_Wtime();

for(int it=0;it<T;it++)
{

if(myrank==0)
{

    for(int i=0;i<pr;i++)
    {
        for(int j=pr-1;j<pr;j++)
        {
            MPI_Isend(&mat[i][j] , 1, MPI_DOUBLE,myrank+1, i, MPI_COMM_WORLD, &reqs[i]);
        }
    }


for(int i=0;i<pr;i++)
{
    MPI_Irecv(&lb[i], 1, MPI_DOUBLE, 0, i, MPI_COMM_WORLD,&reqs[i]);
    
} 
 


for(int j=0;j<pr;j++)
        {
            MPI_Isend(&mat[pr-1][j] , 1, MPI_DOUBLE,myrank+rw, j, MPI_COMM_WORLD, &reqs[j]);
        }
    
    
    for(int i=0;i<pr;i++)
    {
        MPI_Irecv(&ub[i], 1, MPI_DOUBLE, 0, i, MPI_COMM_WORLD,&reqs[i]);
    }
 

}



if(myrank==rw-1)
{

    for(int i=0;i<pr;i++)
    MPI_Isend(&mat[i][0] , 1, MPI_DOUBLE,myrank-1, i, MPI_COMM_WORLD, &reqs[i]);
    
    for(int i=0;i<pr;i++)
    MPI_Irecv(&rb[i], 1, MPI_DOUBLE,rw-1, i, MPI_COMM_WORLD,&reqs[i]);

for(int i=0;i<pr;i++)
    {
       
            MPI_Isend(&mat[pr-1][i] , 1, MPI_DOUBLE,myrank+rw, i, MPI_COMM_WORLD, &reqs[i]);
        
    }


for(int i=0;i<pr;i++)
{
    MPI_Irecv(&ub[i], 1, MPI_DOUBLE, rw-1, i, MPI_COMM_WORLD,&reqs[i]);
    
} 
 

}


if(myrank >0 && myrank<rw-1)
{


    for(int i=0;i<pr;i++)
    MPI_Isend(&mat[i][0] , 1, MPI_DOUBLE,myrank-1, i, MPI_COMM_WORLD, &reqs[i]);
    
    for(int i=0;i<pr;i++)
    MPI_Irecv(&rb[i], 1, MPI_DOUBLE,myrank, i, MPI_COMM_WORLD,&reqs[i]);
 

for(int i=0;i<pr;i++)
    {
        for(int j=pr-1;j<pr;j++)
        {
            MPI_Isend(&mat[i][j] , 1, MPI_DOUBLE,myrank+1, i, MPI_COMM_WORLD, &reqs[i]);
        }
    }


for(int i=0;i<pr;i++)
{
    MPI_Irecv(&lb[i], 1, MPI_DOUBLE,myrank, i, MPI_COMM_WORLD,&reqs[i]);
    
} 
 

for(int i=0;i<pr;i++)
    {
       
            MPI_Isend(&mat[pr-1][i] , 1, MPI_DOUBLE,myrank+rw, i, MPI_COMM_WORLD, &reqs[i]);
        
    }


for(int i=0;i<pr;i++)
{
    MPI_Irecv(&ub[i], 1, MPI_DOUBLE, myrank, i, MPI_COMM_WORLD,&reqs[i]);
    
} 
 

}


if(myrank+rw==p)
{

for(int i=0;i<pr;i++)
    {
        
            MPI_Isend(&mat[0][i] , 1, MPI_DOUBLE,myrank-rw, i, MPI_COMM_WORLD, &reqs[i]);
        
    }
    
    for(int i=0;i<pr;i++)
    MPI_Irecv(&db[i], 1, MPI_DOUBLE,myrank, i, MPI_COMM_WORLD,&reqs[i]);
    
    
for(int i=0;i<pr;i++)
    {
        for(int j=pr-1;j<pr;j++)
        {
            MPI_Isend(&mat[i][j] , 1, MPI_DOUBLE,myrank+1, i, MPI_COMM_WORLD, &reqs[i]);
        }
    }


for(int i=0;i<pr;i++)
{
    MPI_Irecv(&lb[i], 1, MPI_DOUBLE,myrank, i, MPI_COMM_WORLD,&reqs[i]);
    
} 
 

}


if(myrank==p-1)
{

for(int i=0;i<pr;i++)
    MPI_Isend(&mat[i][0] , 1, MPI_DOUBLE,myrank-1, i, MPI_COMM_WORLD, &reqs[i]);
    
    for(int i=0;i<pr;i++)
    MPI_Irecv(&rb[i], 1, MPI_DOUBLE,myrank, i, MPI_COMM_WORLD,&reqs[i]);
 
for(int j=0;j<pr;j++)
        {
            MPI_Isend(&mat[0][j] , 1, MPI_DOUBLE,myrank-rw, j, MPI_COMM_WORLD, &reqs[j]);
        }
    
    
    for(int i=0;i<pr;i++)
    {
        MPI_Irecv(&db[i], 1, MPI_DOUBLE, myrank, i, MPI_COMM_WORLD,&reqs[i]);
    }
 

}


if((myrank+rw)>p && myrank <(p-1))
{


    for(int i=0;i<pr;i++)
    MPI_Isend(&mat[i][0] , 1, MPI_DOUBLE,myrank-1, i, MPI_COMM_WORLD, &reqs[i]);
    
    for(int i=0;i<pr;i++)
    MPI_Irecv(&rb[i], 1, MPI_DOUBLE, myrank, i, MPI_COMM_WORLD,&reqs[i]);
 

for(int i=0;i<pr;i++)
    {
        for(int j=pr-1;j<pr;j++)
        {
            MPI_Isend(&mat[i][j] , 1, MPI_DOUBLE,myrank+1, i, MPI_COMM_WORLD, &reqs[i]);
        }
    }


for(int i=0;i<pr;i++)
{
    MPI_Irecv(&lb[i], 1, MPI_DOUBLE, myrank, i, MPI_COMM_WORLD,&reqs[i]);
    
} 
 
for(int i=0;i<pr;i++)
    {
        
            MPI_Isend(&mat[0][i] , 1, MPI_DOUBLE,myrank-rw, i, MPI_COMM_WORLD, &reqs[i]);
        
    }
    
    for(int i=0;i<pr;i++)
    MPI_Irecv(&db[i], 1, MPI_DOUBLE,myrank, i, MPI_COMM_WORLD,&reqs[i]);
    
    
    
}


if(myrank %rw==0 && myrank!=0 && myrank !=(p-rw))
{

for(int i=0;i<pr;i++)
    {
        
            MPI_Isend(&mat[0][i] , 1, MPI_DOUBLE,myrank-rw, i, MPI_COMM_WORLD, &reqs[i]);
        
    }
    
    for(int i=0;i<pr;i++)
    MPI_Irecv(&db[i], 1, MPI_DOUBLE, myrank, i, MPI_COMM_WORLD,&reqs[i]);
    
    

for(int i=0;i<pr;i++)
    {
        for(int j=pr-1;j<pr;j++)
        {
            MPI_Isend(&mat[i][j] , 1, MPI_DOUBLE,myrank+1, i, MPI_COMM_WORLD, &reqs[i]);
        }
    }


for(int i=0;i<pr;i++)
{
    MPI_Irecv(&lb[i], 1, MPI_DOUBLE,myrank, i, MPI_COMM_WORLD,&reqs[i]);
    
} 
 


        for(int j=0;j<pr;j++)
        {
            MPI_Isend(&mat[pr-1][j] , 1, MPI_DOUBLE,myrank+rw, j, MPI_COMM_WORLD, &reqs[j]);
        }
    
    
    for(int i=0;i<pr;i++)
    {
        MPI_Irecv(&ub[i], 1, MPI_DOUBLE,myrank, i, MPI_COMM_WORLD,&reqs[i]);
    }
 
}


if(myrank%rw==rw-1 && myrank!=rw-1 && myrank!=p-1)
{
 for(int i=0;i<pr;i++)
    MPI_Isend(&mat[i][0] , 1, MPI_DOUBLE,myrank-1, i, MPI_COMM_WORLD, &reqs[i]);
    
    for(int i=0;i<pr;i++)
    MPI_Irecv(&rb[i], 1, MPI_DOUBLE, myrank, i, MPI_COMM_WORLD,&reqs[i]);
 

for(int j=0;j<pr;j++)
        {
            MPI_Isend(&mat[pr-1][j] , 1, MPI_DOUBLE,myrank+rw, j, MPI_COMM_WORLD, &reqs[j]);
        }
    
    
    for(int i=0;i<pr;i++)
    {
        MPI_Irecv(&ub[i], 1, MPI_DOUBLE,myrank, i, MPI_COMM_WORLD,&reqs[i]);
    }


for(int i=0;i<pr;i++)
    MPI_Isend(&mat[0][i] , 1, MPI_DOUBLE,myrank-rw, i, MPI_COMM_WORLD, &reqs[i]);
    
    for(int i=0;i<pr;i++)
    MPI_Irecv(&db[i], 1, MPI_DOUBLE,myrank, i, MPI_COMM_WORLD,&reqs[i]);

}

if(myrank >rw && myrank%rw!=rw-1 && myrank%rw!=0 && myrank<(p-rw))
{
for(int i=0;i<pr;i++)
    {
        
            MPI_Isend(&mat[0][i] , 1, MPI_DOUBLE,myrank-rw, i, MPI_COMM_WORLD, &reqs[i]);
        
    }
    
    for(int i=0;i<pr;i++)
    MPI_Irecv(&db[i], 1, MPI_DOUBLE,myrank, i, MPI_COMM_WORLD,&reqs[i]);
    
    

for(int j=0;j<pr;j++)
        {
            MPI_Isend(&mat[pr-1][j] , 1, MPI_DOUBLE,myrank+rw, j, MPI_COMM_WORLD, &reqs[j]);
        }
    
    
    for(int i=0;i<pr;i++)
    {
        MPI_Irecv(&ub[i], 1, MPI_DOUBLE, myrank, i, MPI_COMM_WORLD,&reqs[i]);
    }


for(int i=0;i<pr;i++)
    MPI_Isend(&mat[i][0] , 1, MPI_DOUBLE,myrank-1, i, MPI_COMM_WORLD, &reqs[i]);
    
    for(int i=0;i<pr;i++)
    MPI_Irecv(&rb[i], 1, MPI_DOUBLE, myrank, i, MPI_COMM_WORLD,&reqs[i]);


for(int i=0;i<pr;i++)
    {
        for(int j=pr-1;j<pr;j++)
        {
            MPI_Isend(&mat[i][j] , 1, MPI_DOUBLE,myrank+1, i, MPI_COMM_WORLD, &reqs[i]);
        }
    }


for(int i=0;i<pr;i++)
{
    MPI_Irecv(&lb[i], 1, MPI_DOUBLE,myrank, i, MPI_COMM_WORLD,&reqs[i]);
    
} 


}




for(int i=1;i<pr-1;i++)
{
    mat[i][0]=0.25*(mat[i-1][0]+mat[i+1][0]+mat[i][1]+lb[i]);
    
}

for(int i=1;i<pr-1;i++)
{
    for(int j=pr-1;j<pr;j++)
    {
        mat[i][j]=0.25*(mat[i-1][j]+mat[i+1][j]+mat[i][j-1]+rb[i]);
    }
}
for(int i=1;i<pr-1;i++)
{
    mat[0][i]=0.25*(mat[1][i]+mat[0][i+1]+ub[i]+mat[0][i-1]);
}

for(int i=pr-1;i<pr;i++)
{
    for(int j=1;j<pr-1;j++)
    {
        mat[i][j]=0.25*(mat[i][j-1]+mat[i][j+1]+mat[i-1][j]+db[j]);
    }
}



for(int i=1;i<pr-1;i++)
{
    for(int j=1;j<pr-1;j++)
    {
        mat[i][j]=0.25*(mat[i][j-1]+mat[i][j+1]+mat[i-1][j]+mat[i+1][j]);
    }
}


mat[0][0]=0.25*(ub[0]+lb[0]+mat[1][0]+mat[0][1]);

mat[0][pr-1]=0.25*(ub[pr-1]+rb[0]+mat[0][pr-2]+mat[1][pr-1]);

mat[pr-1][0]=0.25*(lb[pr-1]+mat[pr-2][0]+mat[pr-1][1]+db[0]);

mat[pr-1][pr-1]=0.25*(db[pr-1]+rb[pr-1]+mat[pr-2][pr-1]+mat[pr-1][pr-2]);

MPI_Barrier(MPI_COMM_WORLD);
}

eTime=MPI_Wtime();
time=eTime-sTime;

if(time>max)
max=time;

if(myrank==0)
printf("%lf\n",max);


free(mat[0]);
free(mat);

MPI_Finalize();

}


if(fun==2)
{
    int myrank, size,N,rw,c,pr,T;
int p;
double *db,*rb,*lb,*ub,*buff;
int position=0;
double max=0.0;

  MPI_Status status;
double sTime,eTime,timee;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  
  p=size;
  N=atoi (argv[1]);
  T=atoi (argv[2]);
  
  rw=sqrt(p);
  pr=sqrt(N);
 
 double **mat;
 
 mat=malloc(pr*sizeof(double*));
 mat[0]=malloc(pr*pr*sizeof(double));
 
 for(int i=1;i<pr;i++)
 mat[i]=mat[0]+pr*i;
 
 
 for(int i=0;i<pr;i++)
 {
     for(int j=0;j<pr;j++)
     {
         mat[i][j]=rand()%10;
     }
 }
 

lb= (double *) malloc(pr * sizeof(double));   
rb=(double *) malloc(pr * sizeof(double));
ub=(double *) malloc(pr * sizeof(double));
db=(double *) malloc(pr * sizeof(double));

 memset(lb, 0, pr * sizeof(double));
 memset(rb, 0, pr * sizeof(double));
 memset(ub, 0, pr * sizeof(double));
 memset(db, 0, pr * sizeof(double));
 
 int prr=N*p*4 *sizeof(double);

buff=(double *) malloc(N*p*4 * sizeof(double));
memset(buff, 0, N*p*4 *sizeof(double));
 
MPI_Request reqs[p];

sTime=MPI_Wtime();

for(int it=0;it<T;it++)
{

if(myrank==0)
{

    for(int i=0;i<pr;i++)
    {
        for(int j=pr-1;j<pr;j++)
        {

MPI_Pack(&mat[i][j] ,1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
            
        }
    }
MPI_Isend(buff, pr, MPI_PACKED,myrank+1, 1, MPI_COMM_WORLD,reqs);


    MPI_Irecv(&lb, pr, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD,reqs);
    



for(int j=0;j<pr;j++)
        {

MPI_Pack(&mat[pr-1][j],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
            
        }
    MPI_Isend(buff, pr, MPI_PACKED,myrank+rw, 1, MPI_COMM_WORLD,reqs);
    
    
        MPI_Irecv(&ub, pr, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD,reqs);
   

}



if(myrank==rw-1)
{

    for(int i=0;i<pr;i++)
MPI_Pack(&mat[i][0],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
    MPI_Isend( buff, pr, MPI_PACKED,myrank-1, 1, MPI_COMM_WORLD,reqs);
    
    
    MPI_Irecv(&rb, pr, MPI_DOUBLE,rw-1, 1, MPI_COMM_WORLD,reqs);
 

for(int i=0;i<pr;i++)
    {
        
MPI_Pack(&mat[pr-1][i],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
            
        
    }
MPI_Isend(buff , pr, MPI_PACKED,myrank+rw, 1, MPI_COMM_WORLD,reqs);


    MPI_Irecv(&ub, pr, MPI_DOUBLE, rw-1, 1, MPI_COMM_WORLD,reqs);
    


}


if(myrank >0 && myrank<rw-1)
{


    for(int i=0;i<pr;i++)
MPI_Pack(&mat[i][0],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
    MPI_Isend(buff , pr, MPI_PACKED,myrank-1, 1, MPI_COMM_WORLD,reqs);
    
    
    MPI_Irecv(&rb, pr, MPI_DOUBLE,myrank, 1, MPI_COMM_WORLD,reqs);
 
 for(int i=0;i<pr;i++)
 {
for(int j=pr-1;j<pr;j++)
    {
        
MPI_Pack(&mat[i][j],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
            
    }    
    }

MPI_Isend(buff , pr, MPI_PACKED,myrank+1, 1, MPI_COMM_WORLD,reqs);

    MPI_Irecv(&lb, pr, MPI_DOUBLE,myrank, 1, MPI_COMM_WORLD,reqs);
    



for(int i=0;i<pr;i++)
    {
       
MPI_Pack(&mat[pr-1][i],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
           
        
    }
 MPI_Isend(buff , pr, MPI_PACKED,myrank+rw, 1, MPI_COMM_WORLD,reqs);


    MPI_Irecv(&ub, pr, MPI_DOUBLE, myrank, 1, MPI_COMM_WORLD,reqs);
    


}


if(myrank+rw==p)
{

for(int i=0;i<pr;i++)
    {
        MPI_Pack(&mat[0][i],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
           
        
    }
     MPI_Isend(buff, pr, MPI_PACKED,myrank-rw, 1, MPI_COMM_WORLD,reqs);
    
    MPI_Irecv(&db, pr, MPI_DOUBLE,myrank, 1, MPI_COMM_WORLD,reqs);
 
for(int i=0;i<pr;i++)
    {
        for(int j=pr-1;j<pr;j++)
        {
MPI_Pack(&mat[i][j],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
           
        }
    }

 MPI_Isend( buff, pr,MPI_PACKED,myrank+1, 1, MPI_COMM_WORLD,reqs);

    MPI_Irecv(&lb, pr, MPI_DOUBLE,myrank, 1, MPI_COMM_WORLD,reqs);
    


}


if(myrank==p-1)
{

for(int i=0;i<pr;i++)
MPI_Pack(&mat[i][0],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
    MPI_Isend( buff, pr, MPI_PACKED,myrank-1, 1, MPI_COMM_WORLD,reqs);
    
    
    MPI_Irecv(&rb, pr, MPI_DOUBLE,myrank, 1, MPI_COMM_WORLD,reqs);
 

    for(int j=0;j<pr;j++)
        {
MPI_Pack(&mat[0][j],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
            
        }
    
    MPI_Isend(buff, pr, MPI_PACKED,myrank-rw, 1, MPI_COMM_WORLD,reqs);
    
        MPI_Irecv(&db, pr, MPI_DOUBLE, myrank, 1, MPI_COMM_WORLD,reqs);
   

}


if((myrank+rw)>p && myrank !=p-1)
{


    for(int i=0;i<pr;i++)
MPI_Pack(&mat[i][0] ,1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
    MPI_Isend(buff, pr,MPI_PACKED,myrank-1, 1, MPI_COMM_WORLD,reqs);
    
    
    MPI_Irecv(&rb, pr, MPI_DOUBLE, myrank, 1, MPI_COMM_WORLD,reqs);

for(int i=0;i<pr;i++)
    {
        for(int j=pr-1;j<pr;j++)
        {
MPI_Pack(&mat[i][j],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
           
        }
    }

 MPI_Isend( buff, pr, MPI_PACKED,myrank+1, 1, MPI_COMM_WORLD,reqs);

    MPI_Irecv(&lb, pr, MPI_DOUBLE, myrank, 1, MPI_COMM_WORLD,reqs);
    


for(int i=0;i<pr;i++)
    {
        MPI_Pack(&mat[0][i],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
    }
        MPI_Isend( buff, pr,MPI_PACKED,myrank-rw, 1, MPI_COMM_WORLD,reqs);
        
    
    
    
    MPI_Irecv(&db, pr, MPI_DOUBLE,myrank, 1, MPI_COMM_WORLD,reqs);
    
   
    
}


if(myrank %rw==0 && myrank!=0 && myrank !=p-rw)
{

for(int i=0;i<pr;i++)
    {
MPI_Pack(&mat[0][i],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
        
            
        
    }
MPI_Isend(buff , pr, MPI_PACKED,myrank-rw, 1, MPI_COMM_WORLD,reqs);
    
    
    MPI_Irecv(&db, pr, MPI_DOUBLE, myrank, 1, MPI_COMM_WORLD,reqs);
    
    
    for(int i=0;i<pr;i++)
    {
        for(int j=pr-1;j<pr;j++)
        {
MPI_Pack(&mat[i][j],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
           
        }
    }
    

 MPI_Isend(buff , pr, MPI_PACKED,myrank+1, 1, MPI_COMM_WORLD,reqs);

    MPI_Irecv(&lb, pr, MPI_DOUBLE,myrank, 1, MPI_COMM_WORLD,reqs);
    



        for(int j=0;j<pr;j++)
        {
MPI_Pack(&mat[pr-1][j],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
           
        }
    
     MPI_Isend(buff,pr, MPI_PACKED,myrank+rw, 1, MPI_COMM_WORLD,reqs);
    
        MPI_Irecv(&ub, pr, MPI_DOUBLE,myrank, 1, MPI_COMM_WORLD,reqs);
    

}



if(myrank >rw && myrank%rw!=rw-1 && myrank%rw!=0 && myrank<(p-rw))
{
for(int i=0;i<pr;i++)
    {
MPI_Pack(&mat[0][i],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
        
            
        
    }
    
MPI_Isend(buff ,pr, MPI_PACKED,myrank-rw, 1, MPI_COMM_WORLD,reqs);
   
    MPI_Irecv(&db, pr, MPI_DOUBLE,myrank, 1, MPI_COMM_WORLD,reqs);
    
    
for(int j=0;j<pr;j++)
        {
MPI_Pack(&mat[pr-1][j],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
          }
  MPI_Isend(buff,pr, MPI_PACKED,myrank+rw, 1, MPI_COMM_WORLD,reqs);
        
    
    
   
        MPI_Irecv(&ub, pr, MPI_DOUBLE, myrank, 1, MPI_COMM_WORLD,reqs);
    

for(int i=0;i<pr;i++)
MPI_Pack(&mat[i][0],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
    MPI_Isend(buff ,pr, MPI_PACKED,myrank-1,1, MPI_COMM_WORLD,reqs);
    
    
    MPI_Irecv(&rb,pr, MPI_DOUBLE, myrank, 1, MPI_COMM_WORLD,reqs);
 

for(int i=0;i<pr;i++)
    {
        for(int j=pr-1;j<pr;j++)
        {
MPI_Pack(&mat[i][j],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
            
        }
    }

MPI_Isend( buff,pr, MPI_PACKED,myrank+1, 1, MPI_COMM_WORLD,reqs);

    MPI_Irecv(&lb, pr, MPI_DOUBLE,myrank, 1, MPI_COMM_WORLD,reqs);
    
}


if(myrank%rw==rw-1 && myrank!=rw-1 && myrank!=p-1)
{
    
    for(int i=0;i<pr;i++)
    {
MPI_Pack(&mat[0][i],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
        
            
        
    }
    
MPI_Isend(buff ,pr, MPI_PACKED,myrank-rw, 1, MPI_COMM_WORLD,reqs);
   
    MPI_Irecv(&db, pr, MPI_DOUBLE,myrank, 1, MPI_COMM_WORLD,reqs);
    
    
    for(int i=0;i<pr;i++)
    {
       
MPI_Pack(&mat[pr-1][i],1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
           
        
    }
 MPI_Isend(buff , pr, MPI_PACKED,myrank+rw, 1, MPI_COMM_WORLD,reqs);


    MPI_Irecv(&ub, pr, MPI_DOUBLE, myrank, 1, MPI_COMM_WORLD,reqs);
    

for(int i=0;i<pr;i++)
MPI_Pack(&mat[i][0] ,1,MPI_DOUBLE,buff,prr,&position,MPI_COMM_WORLD);
    MPI_Isend(buff, pr,MPI_PACKED,myrank-1, 1, MPI_COMM_WORLD,reqs);
    
    
    MPI_Irecv(&rb, pr, MPI_DOUBLE, myrank, 1, MPI_COMM_WORLD,reqs);

    
    
}

for(int i=1;i<pr-1;i++)
{
    mat[i][0]=0.25*(mat[i-1][0]+mat[i+1][0]+mat[i][1]+lb[i]);
    
}

for(int i=1;i<pr-1;i++)
{
    for(int j=pr-1;j<pr;j++)
    {
        mat[i][j]=0.25*(mat[i-1][j]+mat[i+1][j]+mat[i][j-1]+rb[i]);
    }
}
for(int i=1;i<pr-1;i++)
{
    mat[0][i]=0.25*(mat[1][i]+mat[0][i+1]+ub[i]+mat[0][i-1]);
}

for(int i=pr-1;i<pr;i++)
{
    for(int j=1;j<pr-1;j++)
    {
        mat[i][j]=0.25*(mat[i][j-1]+mat[i][j+1]+mat[i-1][j]+db[j]);
    }
}



for(int i=1;i<pr-1;i++)
{
    for(int j=1;j<pr-1;j++)
    {
        mat[i][j]=0.25*(mat[i][j-1]+mat[i][j+1]+mat[i-1][j]+mat[i+1][j]);
    }
}


mat[0][0]=0.25*(ub[0]+lb[0]+mat[1][0]+mat[0][1]);

mat[0][pr-1]=0.25*(ub[pr-1]+rb[0]+mat[0][pr-2]+mat[1][pr-1]);

mat[pr-1][0]=0.25*(lb[pr-1]+mat[pr-2][0]+mat[pr-1][1]+db[0]);

mat[pr-1][pr-1]=0.25*(db[pr-1]+rb[pr-1]+mat[pr-2][pr-1]+mat[pr-1][pr-2]);

MPI_Barrier(MPI_COMM_WORLD);
}


eTime=MPI_Wtime();
timee=eTime-sTime;

if(timee>max)
max=timee;

if(myrank==0)
printf("%lf\n",max);

free(mat[0]);
free(mat);

MPI_Finalize();
}

if(fun==3)
{
    int myrank, size,N,rw,c,pr,T;
int p;
double *db,*rb,*lb,*ub;

  MPI_Status status;
double sTime,eTime,time,max=0.0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  
  p=size;
  N=atoi (argv[1]);
  T=atoi (argv[2]);
  
  rw=sqrt(p);
  pr=sqrt(N);
  
  
  MPI_Datatype ctype,vtype;
  
  MPI_Type_contiguous(pr,MPI_DOUBLE,&ctype);
  MPI_Type_commit(&ctype);
  
  MPI_Type_vector(pr,1,pr,MPI_DOUBLE,&vtype);
   MPI_Type_commit(&vtype);
 
double **mat;
 
 mat=malloc(pr*sizeof(double*));
 mat[0]=malloc(pr*pr*sizeof(double));
 
 for(int i=1;i<pr;i++)
 mat[i]=mat[0]+pr*i;
 
 
 for(int i=0;i<pr;i++)
 {
     for(int j=0;j<pr;j++)
     {
         mat[i][j]=rand()%10;
     }
 }
 

lb= (double *) malloc(pr * sizeof(double));   
rb=(double *) malloc(pr * sizeof(double));
ub=(double *) malloc(pr * sizeof(double));
db=(double *) malloc(pr * sizeof(double));

 memset(lb, 0, pr * sizeof(double));
 memset(rb, 0, pr * sizeof(double));
 memset(ub, 0, pr * sizeof(double));
 memset(db, 0, pr * sizeof(double));

 
MPI_Request reqs[p];

sTime=MPI_Wtime();

for(int it=0;it<T;it++)
{

if(myrank==0)
{

    
            MPI_Isend(&mat[0][pr-1] , 1, vtype,myrank+1, 1, MPI_COMM_WORLD, reqs);
        
    



    MPI_Irecv(&lb, pr, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD,reqs);
    

 


for(int j=0;j<pr;j++)
        {
            MPI_Isend(&mat[pr-1][0] , 1, ctype,myrank+rw, 2, MPI_COMM_WORLD, reqs);
        }
    
    
   
        MPI_Irecv(&ub, pr, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD,reqs);
    
 

}



if(myrank==rw-1)
{

    
    MPI_Isend(&mat[0][0] , 1, vtype,myrank-1, 3, MPI_COMM_WORLD, reqs);
    
   
    MPI_Irecv(&rb, pr, MPI_DOUBLE,rw-1, 3, MPI_COMM_WORLD,reqs);


       
            MPI_Isend(&mat[pr-1][0] , 1, ctype,myrank+rw, 4, MPI_COMM_WORLD, reqs);
        
    



    MPI_Irecv(&ub, pr, MPI_DOUBLE, rw-1, 4, MPI_COMM_WORLD,reqs);
    

 

}


if(myrank >0 && myrank<rw-1)
{


    
    MPI_Isend(&mat[0][0] , 1, vtype,myrank-1, 5, MPI_COMM_WORLD, reqs);
    
    
    MPI_Irecv(&rb, pr, MPI_DOUBLE,myrank, 5, MPI_COMM_WORLD,reqs);
 


            MPI_Isend(&mat[0][pr-1] , 1, vtype,myrank+1, 6, MPI_COMM_WORLD, reqs);
       


    MPI_Irecv(&lb, pr, MPI_DOUBLE,myrank, 6, MPI_COMM_WORLD,reqs);
    

 


       
            MPI_Isend(&mat[pr-1][0] , 1, ctype,myrank+rw, 7, MPI_COMM_WORLD, reqs);
        
    



    MPI_Irecv(&ub, pr, MPI_DOUBLE, myrank, 7, MPI_COMM_WORLD,reqs);
    

 

}


if(myrank+rw==p)
{


        
            MPI_Isend(&mat[0][0] , 1, ctype,myrank-rw, 8, MPI_COMM_WORLD, reqs);
        
    
    
    
    MPI_Irecv(&db, pr, MPI_DOUBLE,myrank, 8, MPI_COMM_WORLD,reqs);
    
    

            MPI_Isend(&mat[0][pr-1] , 1, vtype,myrank+1, 9, MPI_COMM_WORLD, reqs);
        


    MPI_Irecv(&lb, pr, MPI_DOUBLE,myrank, 9, MPI_COMM_WORLD,reqs);
    

 

}


if(myrank==p-1)
{


    MPI_Isend(&mat[0][0] , 1, vtype,myrank-1, 10, MPI_COMM_WORLD, reqs);
    
    
    MPI_Irecv(&rb, pr, MPI_DOUBLE,myrank, 10, MPI_COMM_WORLD,reqs);
 

            MPI_Isend(&mat[0][pr-1] , 1, vtype,myrank-rw, 11, MPI_COMM_WORLD, reqs);
        
    
    
   
        MPI_Irecv(&db, pr, MPI_DOUBLE, myrank, 11, MPI_COMM_WORLD,reqs);
   
 

}


if((myrank+rw)>p && myrank !=(p-1))
{


    
    MPI_Isend(&mat[0][0] , 1, vtype,myrank-1, 12, MPI_COMM_WORLD, reqs);
    
    
    MPI_Irecv(&rb, pr, MPI_DOUBLE, myrank, 12, MPI_COMM_WORLD,reqs);
 


            MPI_Isend(&mat[0][pr-1] , 1, vtype,myrank+1, 13, MPI_COMM_WORLD, reqs);
       



    MPI_Irecv(&lb, pr, MPI_DOUBLE, myrank, 13, MPI_COMM_WORLD,reqs);
    

 

        
            MPI_Isend(&mat[0][0] , 1, ctype,myrank-rw, 14, MPI_COMM_WORLD, reqs);
        
    
    
    
    MPI_Irecv(&db, pr, MPI_DOUBLE,myrank, 14, MPI_COMM_WORLD,reqs);
    
    
    
}


if(myrank %rw==0 && myrank!=0 && myrank !=(p-rw))
{

for(int i=0;i<pr;i++)
    {
        
            MPI_Isend(&mat[0][0] , 1, ctype,myrank-rw, 15, MPI_COMM_WORLD, reqs);
        
    }
    
    
    MPI_Irecv(&db, pr, MPI_DOUBLE, myrank, 15, MPI_COMM_WORLD,reqs);
    
    


            MPI_Isend(&mat[0][pr-1] , 1, vtype,myrank+1, 16, MPI_COMM_WORLD, reqs);
        



    MPI_Irecv(&lb, pr, MPI_DOUBLE,myrank, 16, MPI_COMM_WORLD,reqs);
    

 


        
            MPI_Isend(&mat[pr-1][0] , 1, ctype,myrank+rw, 18, MPI_COMM_WORLD, reqs);
        
    
    
    
        MPI_Irecv(&ub, pr, MPI_DOUBLE,myrank, 18, MPI_COMM_WORLD,reqs);
    
 
}


if(myrank%rw==rw-1 && myrank!=rw-1 && myrank!=p-1)
{
 
    MPI_Isend(&mat[0][0] , 1, vtype,myrank-1, 19, MPI_COMM_WORLD, reqs);
    
    
    MPI_Irecv(&rb, pr, MPI_DOUBLE, myrank, 19, MPI_COMM_WORLD,reqs);
 


            MPI_Isend(&mat[pr-1][0] , 1, ctype,myrank+rw, 19, MPI_COMM_WORLD, reqs);
        
    
    
    
        MPI_Irecv(&ub, pr, MPI_DOUBLE,myrank, 19, MPI_COMM_WORLD,reqs);
    



    MPI_Isend(&mat[0][0] , 1, ctype,myrank-rw, 20, MPI_COMM_WORLD, reqs);
    
    
    MPI_Irecv(&db, pr, MPI_DOUBLE,myrank, 20, MPI_COMM_WORLD,reqs);

}

if(myrank >rw && myrank%rw!=rw-1 && myrank%rw!=0 && myrank<(p-rw))
{

        
            MPI_Isend(&mat[0][0] , 1, ctype,myrank-rw, 21, MPI_COMM_WORLD, reqs);
        
    
    
    
    MPI_Irecv(&db, pr, MPI_DOUBLE,myrank, 21, MPI_COMM_WORLD,reqs);
    
    


            MPI_Isend(&mat[pr-1][0] , 1, ctype,myrank+rw, 22, MPI_COMM_WORLD, reqs);
        
    
    
    
        MPI_Irecv(&ub, pr, MPI_DOUBLE, myrank, 22, MPI_COMM_WORLD,reqs);
    



    MPI_Isend(&mat[0][0] , 1, vtype,myrank-1, 23, MPI_COMM_WORLD, reqs);
    
    
    MPI_Irecv(&rb, pr, MPI_DOUBLE, myrank, 23, MPI_COMM_WORLD,reqs);



            MPI_Isend(&mat[0][pr-1] , 1, vtype,myrank+1, 24, MPI_COMM_WORLD, reqs);
        



    MPI_Irecv(&lb, pr, MPI_DOUBLE,myrank, 24, MPI_COMM_WORLD,reqs);
    
 


}



for(int i=1;i<pr-1;i++)
{
    mat[i][0]=0.25*(mat[i-1][0]+mat[i+1][0]+mat[i][1]+lb[i]);
    
}

for(int i=1;i<pr-1;i++)
{
    for(int j=pr-1;j<pr;j++)
    {
        mat[i][j]=0.25*(mat[i-1][j]+mat[i+1][j]+mat[i][j-1]+rb[i]);
    }
}
for(int i=1;i<pr-1;i++)
{
    mat[0][i]=0.25*(mat[1][i]+mat[0][i+1]+ub[i]+mat[0][i-1]);
}

for(int i=pr-1;i<pr;i++)
{
    for(int j=1;j<pr-1;j++)
    {
        mat[i][j]=0.25*(mat[i][j-1]+mat[i][j+1]+mat[i-1][j]+db[j]);
    }
}



for(int i=1;i<pr-1;i++)
{
    for(int j=1;j<pr-1;j++)
    {
        mat[i][j]=0.25*(mat[i][j-1]+mat[i][j+1]+mat[i-1][j]+mat[i+1][j]);
    }
}


mat[0][0]=0.25*(ub[0]+lb[0]+mat[1][0]+mat[0][1]);

mat[0][pr-1]=0.25*(ub[pr-1]+rb[0]+mat[0][pr-2]+mat[1][pr-1]);

mat[pr-1][0]=0.25*(lb[pr-1]+mat[pr-2][0]+mat[pr-1][1]+db[0]);

mat[pr-1][pr-1]=0.25*(db[pr-1]+rb[pr-1]+mat[pr-2][pr-1]+mat[pr-1][pr-2]);

MPI_Barrier(MPI_COMM_WORLD);
}

eTime=MPI_Wtime();
time=eTime-sTime;

if(time>max)
max=time;

if(myrank==0)
printf("%lf\n",max);

free(mat[0]);
free(mat);

MPI_Finalize();
}
  return 0;

}
