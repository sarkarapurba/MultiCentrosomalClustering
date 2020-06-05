/*******************************************************************************/
/*    Mechanics of multi-centrosomal clustering in bipolar mitotic spindles    */
/*									       */									
/*			THIS CODE USES                   		       */
/*      	     MONTE CARLO ALGORITHM           			       */ 
/*      	   TO GENERATE THE MECHANICAL        			       */
/*      	  EQUILIBRIUM CONFIGURATION OF     			       */
/*      	    SPINDLE DURING MITOSIS  				       */
/*******************************************************************************/
/*			DESIGN & DEVELOPMENTS 				       */
/*									       */
/*			S Chatterjee, A Sarkar, R Paul                         */
/*	Indian Association for the Cultivation of Science, Kolkata, India      */
/*      			                                               */
/*				A Mogilner                                     */
/*   	Courant Institute, New York University, New York, USA		       */
/*   				(c) 2020                                       */
/*                                                                             */
/* *****************************************************************************/         



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/**********dyanamic memory allocation********/  
double   *vector ( int nrl, int nrh);
double  **matrix ( int nrl, int nrh, int ncl , int nch );
double ***tensor3 ( int nxl, int nxh, int nyl, int nyh, int nzl, int nzh );

int   *ivector ( int nrl, int nrh );
int  **imatrix ( int nrl, int nrh, int ncl, int nch );
int ***itensor3 ( int nxl, int nxh, int nyl, int nyh, int nzl, int nzh );


double *vector( int nl, int nh)
{
        double *v;
	
        v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
        if (!v) printf("allocation failure in vector()\n\n");
        return v-nl;
}

double **matrix ( int nrl, int nrh, int ncl, int nch )
{
  int i;
  double **m;

  m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
  if (!m) printf("allocation failure 1 in matrix()\n\n");
  m -= nrl;
  
  for(i=nrl;i<=nrh;i++) {
    m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
    if (!m[i]) printf("allocation failure 2 in matrix()\n\n");
    m[i] -= ncl;
  }
  return m;
}

double ***tensor3 ( int nxl, int nxh, int nyl, int nyh, int nzl, int nzh )
{
  int i, j;
  double ***m;
  
  m=(double ***) malloc((unsigned) (nxh-nxl+1)*sizeof(double*));
  if (!m) printf("allocation failure 1 in tensor3()\n\n");
  m -= nxl;
  
  for(i=nxl;i<=nxh;i++) {
    m[i]=(double **) malloc((unsigned) (nyh-nyl+1)*sizeof(double*));
    if (!m[i]) printf("allocation failure 2 in tensor3()\n\n");
    m[i] -= nyl;
  };
  
  for(i=nxl;i<=nxh;i++) {
    for(j=nyl;j<=nyh;j++) {
      m[i][j]=(double *) malloc((unsigned) (nzh-nzl+1)*sizeof(double));
      if (!m[i][j]) printf("allocation failure 3 in tensor3()\n\n");
      m[i][j] -= nzl;
    }
  };
  
  return m;
  
}

int *ivector( int nl, int nh)
{
  int *v;
  
  v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
  if (!v) printf("allocation failure in ivector()\n\n");
  return v-nl;
}

int **imatrix ( int nrl, int nrh, int ncl, int nch )
{
  int i;
  int **m;
  
  m=(int **) malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
  if (!m) printf("allocation failure 1 in matrix()\n\n");
  m -= nrl;
  
  for(i=nrl;i<=nrh;i++) {
    m[i]=(int *) malloc((unsigned) (nch-ncl+1)*sizeof(int));
    if (!m[i]) printf("allocation failure 2 in matrix()\n\n");
    m[i] -= ncl;
  }
  return m;
}

int ***itensor3 ( int nxl, int nxh, int nyl, int nyh, int nzl, int nzh )
{
  int i, j;
  int ***m;
  
  m=(int ***) malloc((unsigned) (nxh-nxl+1)*sizeof(int*));
  if (!m) printf("allocation failure 1 in tensor3()\n\n");
  m -= nxl;
  
  for(i=nxl;i<=nxh;i++) {
    m[i]=(int **) malloc((unsigned) (nyh-nyl+1)*sizeof(int*));
    if (!m[i]) printf("allocation failure 2 in tensor3()\n\n");
    m[i] -= nyl;
  };
  
  for(i=nxl;i<=nxh;i++) {
    for(j=nyl;j<=nyh;j++) {
      m[i][j]=(int *) malloc((unsigned) (nzh-nzl+1)*sizeof(int));
      if (!m[i][j]) printf("allocation failure 3 in tensor3()\n\n");
      m[i][j] -= nzl;
    }
  };
  
  return m;
  
}
/************************************************************/


/************************************************************/
/*       RANDOM NUMBER GENERATOR AND                        */
/*        LIST OF SUBROUTINES                               */
/************************************************************/
#define N_BITS 31
#define IRND1  250
#define IRND2  103
#define MAX_RN 0X7FFFFFFF
#define MAX_R2 0X40000000

int ir[256], irndx1[256], irndx2[256];

void rndini (int seed);
double ran2(int *idum);
int Pw(double n);
int OutOfSpheroid_ini(int x, int y, int z);
int OutOfSpheroid(int x, int y, int z);
void DefineCellVol();
void DefineSurface();
void OccConfig();
void NearestNeighbAllocation();
void Monte_Carlo_Update();
void Particle_Exchange(int S, int X_S, int Y_S, int Z_S, int NB, int X_NB, int Y_NB, int Z_NB);
void PrintPNGSnpshot(int n);
void Visualization_povray(int tt, int ensemble);
void Final_configuration(int ensemble,int MAX_POWER);
void IdentifyPoles();
void OutputPolesHist(int ensemble);
void OutputPoles();


///declaration of variables
int idum = -12347;
int counter = 0, t, Lx = 0, Ly, Lz, Lxy, N, MAX_POWER = 0, t1,t2, Nclust,N_total,N_cent,N_ch; 
int MAX_ENSEMBLE = 1, MIN_ENSEMBLE=1, Npoles=0,Npole_inst=0;
int Ax, Ay, Az, X0, Y0, Z0, Ax_ini, Ay_ini, Az_ini;
int *particle, *occ, *label_particle, *particle_dummy;
int  **neighb;
double Temp = 0.5, Lav = 10.0, Lav_ast = 5.0, fcs_csd=0.0, fcs_csk=0.0, fcs_kt=0.0, fcs_ch=0.0, fch_ch=0.0, fkt_kt=0.0, 
  r_cut, fch_wal=0.0, fcs_wal=0.0, fcs_surf=0.0, dist_merge=1.0;



/*******************************************************/  

void nrerror0(const char error_text[])
{
  printf("Numerical Recipes run-time error...\n");
  printf("%s\n",error_text);
  printf("...now exiting to system...\n");
  exit(1);
}
/**********************************************************/ 


/***** random number generator ran2 ***********************/
#include <math.h>

#define M 714025
#define IA 1366
#define IC 150889

double ran2(int *idum)
  
{
  static long iy,irn[98];
  static int iff=0;
  int j;
  void nrerror();
  
  if (*idum < 0 || iff == 0) {
    iff=1;
    if ((*idum=(IC-(*idum)) % M) < 0) *idum = -(*idum);
    for (j=1;j<=97;j++) {
      *idum=(IA*(*idum)+IC) % M;
      irn[j]=(*idum);
    }
    *idum=(IA*(*idum)+IC) % M;
    iy=(*idum);
  }
  j=1 +(int)( 97.0*iy/M);
  if (j > 97 || j < 1) nrerror0("RAN2: This cannot happen.");
  iy=irn[j];
  *idum=(IA*(*idum)+IC) % M;
  irn[j]=(*idum);
  return (double)iy/M;
}

#undef M
#undef IA
#undef IC
/****************************************************************/


/************* 2^n**********************************/ 
int Pw(double n) {
  return((int)(pow(2.0,n)));
}  

///Subroutine for intial placement of the centrosomes and chromosomes within the spheroidal volume having sami-axes: Ax_ini, Ay_ini, Az_ini
int OutOfSpheroid_ini(int x, int y, int z){
  //Here, we assume the spheroid is centered at X0,Y0,Z0 and
  //Ax_ini,Ay_ini,Az_ini are the radii in 3 directions
  
  
  x = x-X0;
  y = y-Y0;
  z = z-Z0;
  
  double x2 = (double)(x*x);
  double y2 = (double)(y*y);
  double z2 = (double)(z*z);
  double Ax2 = (double)(Ax_ini*Ax_ini);
  double Ay2 = (double)(Ay_ini*Ay_ini);
  double Az2 = (double)(Az_ini*Az_ini);
  
  
  if(x2/Ax2+y2/Ay2+z2/Az2<1.0){
    return(0);//inside the spheroid
  }else return(1);//outside the spheroid
}

///Subroutines to confine the particle's movement within the spheroidal volume having semi-axes: Ax, Ay, Az
int OutOfSpheroid(int x, int y, int z){
  //Here, we assume the spheroid is centered at X0,Y0,Z0 and
  //Ax,Ay,Az are the radii in 3 directions
  
  x = x-X0;
  y = y-Y0;
  z = z-Z0;
  
  double x2 = (double)(x*x);
  double y2 = (double)(y*y);
  double z2 = (double)(z*z);
  double Ax2 = (double)(Ax*Ax);
  double Ay2 = (double)(Ay*Ay);
  double Az2 = (double)(Az*Az);
  
  
  if(x2/Ax2+y2/Ay2+z2/Az2<1.0){
    return(0);//inside the spheroid
  }else return(1);//outside the spheroid
}


/// subroutine to define the lattice sites inside the cellular volume 
void DefineCellVol(){
  int i,j,k,site;
  for(k=1;k<=Lz;k++){
    for(i=1;i<=Lx;i++){
      for(j=1;j<=Ly;j++){
	site = j + (i-1)*Lx + (k-1)*Lxy;
	occ[site] = 0;
	particle[site] = 0;
	if(OutOfSpheroid(i,j,k)==1){
	  occ[site] = 2;
	};
      };
    };
  }; 
}

/// Subroutine to define the sufrace nodes of the spheroidal cell///
void DefineSurface(){
  int i,j,k,l,site,nb;
  for(k=1;k<=Lz;k++){
    for(i=1;i<=Lx;i++){
      for(j=1;j<=Ly;j++){
	site = j + (i-1)*Lx + (k-1)*Lxy;
	if(occ[site] == 2){
	  for(l=1;l<=6;l++){
	    nb=neighb[site][l];
	    if(occ[nb] !=2)particle[site]=3;
	  };
	};
      };
    };
  };
}

/**************Particle (centrosomes/chromosomes) Occupancy in the lattice*****************/   
void OccConfig(){
  int i,count;
  //Define Cell Volume
  DefineCellVol();

  //Define Cell Surface
  DefineSurface();
  
  
  
  // distribute all the centrosomes and chromosomes within the volume having semi-axes Ax_ini, Ay_ini, Az_ini
  count = 0;
  if(N_cent>0){
    do{
      int x = X0 + 2*(int)((double)Ax_ini*(0.5-ran2(&idum)));  
      int y = Y0 + 2*(int)((double)Ay_ini*(0.5-ran2(&idum)));
      int z = Z0 + 2*(int)((double)Az_ini*(0.5-ran2(&idum)));
      
      if(x==0)x=1;
      if(y==0)y=1;
      if(z==0)z=1;
      
      if(OutOfSpheroid_ini(x,y,z)==0){
	
	int site = y+(x-1)*Lx+(z-1)*Lxy;
	
	if(occ[site]==0){
	  occ[site]=1;
	  particle[site]= 1;
	  count++;
	};
      };
    }while(count<N_cent);
  };
  
  count = 0;
  if(N_ch>0){
    do{
      int x = X0 + 2*(int)((double)Ax_ini*(0.5-ran2(&idum)));
      int y = Y0 + 2*(int)((double)Ay_ini*(0.5-ran2(&idum)));
      int z = Z0 + 2*(int)((double)Az_ini*(0.5-ran2(&idum)));
      if(x==0)x=1;
      if(y==0)y=1;
      if(z==0)z=1;
      if(OutOfSpheroid_ini(x,y,z)==0){
	int site = y+(x-1)*Lx+(z-1)*Lxy;
	if(occ[site]==0){
	  occ[site]=1;
	  particle[site]= 2;
	  count++;
	};
      };
    }while(count<N_ch);
  };
}

/***********neighbouring array with pbc******************/
void NearestNeighbAllocation(){
  int i,j;
  
  /* neigbor array*/
  for( i=1; i<=N; i++ ){  /*                             4  5  */
    neighb[i][1]=i+1;     /* right      neighbor         |/    */
    neighb[i][2]=i+Lx;     /* bottom     neighbor      3-i-1   */
    neighb[i][3]=i-1;     /* left       neighbor        /|     */
    neighb[i][4]=i-Lx;     /* top        neighbor     6  2     */
    neighb[i][5]=i-Lxy;   /* previous   neighbor               */
    neighb[i][6]=i+Lxy;   /* next       neighbor               */     
  };
  
  /* boundary condition */
  for(i=1;i<=Lx;i++){
    for ( j=1; j<=Ly; j++ ){
      neighb [j+(Lxy*i-Lx)]        [2] = j + (i-1)*Lxy;
      neighb [j+(i-1)*Lxy]         [4] = j + (Lxy*i -Lx);
      
      neighb [j*Lx+(i-1)*Lxy]      [1] = 1+(j-1)*Lx+(i-1)*Lxy;     
      neighb [1+(j-1)*Lx+(i-1)*Lxy] [3] = j*Lx+(i-1)*Lxy;
      
      neighb [j + (i-1)*Lz]         [5] = j + (i-1)*Lz  + N - Lxy;
      neighb [j + (i-1)*Lz + N -Lxy][6] = j + (i-1)*Lz;
    };
  };
}

/**************************************************************/
/******* Monte Carlo Update of particle's position *******/  
void Monte_Carlo_Update(){
  int i,j,k,l,m,n,x,y,z,nb,s;
  int excess;
  double flag;
  
  /* main monte carlo routine */
  for(k=1;k<=Lz;k++){
    for(i=1;i<=Lx;i++){
      for(j=1;j<=Ly;j++){
	s = j + (i-1)*Lx + (k-1)*Lxy;
	x = (int)(double)(s%Lxy)/((double)Lx+0.1)+1;
	if (s%Lxy==0) x=Lx;
	y = (s-1)%Lx+1;
	z = (int)((double)s/((double)Lxy+0.1))+1;
	if(occ[s] == 1){
	  flag = ran2(&idum);
	  nb = (int)(6.0*flag)+1;
	  n = neighb[s][nb];
	  x = (int)(double)(n%Lxy)/((double)Lx+0.1)+1;
	  if (s%Lxy==0) x=Lx;
	  y = (n-1)%Lx+1;
	  z = (int)((double)n/((double)Lxy+0.1))+1;
	  if(occ[n] == 0) {   // Hopping of centrosomes/ chromosomes are allowed if the neighboring site is vacant
	    Particle_Exchange(s,i,j,k,n,x,y,z); 
	  }; 
	};
      };
    };
  };
}


/**************************************************************/
/***Subrountine for the movement of centrosomes and chromosomes to it's randomly chosen neighbor***/
void Particle_Exchange(int S, int X_S, int Y_S, int Z_S, int NB, int X_NB, int Y_NB, int Z_NB){
  
  int i,j,k,l,sum,site;
  double DeltaE,flag,energy1,energy2,prob;
  int sumleft = 0;
  int nb_S,nb_NB;

  DeltaE = 0.0;
  for(k=1;k<=Lz;k++){
    for(i=1;i<=Lx;i++){
      for(j=1;j<=Ly;j++){
	site = j + (i-1)*Lx + (k-1)*Lxy;
	
	if((site != S) && (occ[site]!=0)){
	  double r1 = sqrt((double)((X_S-i)*(X_S-i)+ (Y_S-j)*(Y_S-j) + (Z_S-k)*(Z_S-k))); 
	  double r2 = sqrt((double)((X_NB-i)*(X_NB-i)+ (Y_NB-j)*(Y_NB-j) + (Z_NB-k)*(Z_NB-k)));
	  if(particle[S]==1){//particle[S]=1 is centrosome 
	    if(particle[site]==1){//particle[site]=1 is centrosome
	      DeltaE += fcs_csk*Lav*(((r2*exp(-r2/Lav))-(r1*exp(-r1/Lav)))+(Lav*(exp(-r2/Lav)-exp(-r1/Lav))));
	      DeltaE += fcs_csd*Lav*(exp(-r2/Lav)-exp(-r1/Lav));

	    };
	    if(particle[site]==2){//particle[S]=2 is kinetochore and chromosome
	      
	      DeltaE += -fcs_kt*(r2-r1);
	      DeltaE += fcs_ch*Lav*(exp(-r2/Lav)-exp(-r1/Lav));
	    
	    };
	    
	    if(particle[site]==3){//particle[S]=3 is the boundary, interaction with centrosome
	      DeltaE += fcs_surf*Lav_ast*(exp(-r2/ Lav_ast)-exp(-r1/ Lav_ast));
	      
	    };
	  };
	  
	  if(particle[S]==2){//particle[S]=2 is kinetochore and chromosome
	    if(particle[site]==1){//particle[site]=1 is centrosome
	      DeltaE += -fcs_kt*(r2-r1);
	      DeltaE += fcs_ch*Lav*(exp(-r2/Lav)-exp(-r1/Lav));
	    };
	    
	    //centrosome-boundary steric repulsion
	    energy1=energy2=0.0;
	    for(l=1;l<=6;l++){
	      nb_S = neighb[S][l];
	      if(occ[nb_S]==2)
		energy1 += fcs_wal;
	      nb_NB = neighb[NB][l];
	      if(occ[nb_NB]==2)
		energy2 += fcs_wal;
	    };
	    DeltaE += energy2-energy1;
	      
	  };
	  
	  
	  //chromosome-chromosome steric repulsion
	  if((particle[S]==2) && (particle[site]==2)){//particle[S]=2 is kinetochore and chromosome
	    energy1=energy2=0.0;
	    if(r1<=r_cut) 
	      energy1 = fch_ch*(1.0/r1);
	    if(r2<=r_cut) 
	      energy2 = fch_ch*(1.0/r2);
	    DeltaE += energy2-energy1;
	  };
	  
	  //chromosome-boundary steric repulsion
	  if(particle[S]==2){
	    energy1=energy2=0.0;
	    for(l=1;l<=6;l++){
	      nb_S = neighb[S][l];
	      if(occ[nb_S]==2)
		energy1 += fch_wal;
	      nb_NB = neighb[NB][l];
	      if(occ[nb_NB]==2)
		energy2 += fch_wal;
	    };
	    DeltaE += energy2-energy1;
	  };
	  
	  //kt-kt long distance attraction
	  if((particle[S]==2) && (particle[site]==2)){//particle[S]=2 is kinetochore and chromosome
	    if((r1>4.*r_cut) && (r2>4.*r_cut))
	      DeltaE += fkt_kt*Lav*(exp(-r2/Lav)-exp(-r1/Lav));
	  };
	};
      };
    };
  };
  
  if(DeltaE<=0.0){
    particle[NB] = particle[S];
    particle[S] = 0;
    occ[NB]=1;
    occ[S]=0;
  }else{ 
    flag = ran2(&idum);
    prob = exp(-DeltaE/Temp) ;
    if (flag  < prob) {
      particle[NB] = particle[S];
      particle[S] = 0;
      occ[NB]=1;
      occ[S]=0;
    };
  };
}


/***************** particle configuration ***********/
/************** Output stored in Files ****************/
void PrintPNGSnpshot(int n){
  
  char fn1[100];
  FILE *fp1;
   
  int x,r,k,i,j,time,s;
  double val,denom;
  
  sprintf(fn1,"%dt_%dL_%1.2fT_ChCentSnpshot",1000+n,Lx,Temp);
  fp1 = fopen(fn1,"w");
  fprintf(fp1,"P3\n");
  fprintf(fp1,"#Snapshot of CH-Cent-KT\n");
  fprintf(fp1,"%d  %d\n",Lx,Ly);
  fprintf(fp1,"15\n");
  
  for(k=Lz/2;k<=Lz/2;k++){
    for(i=1;i<=Lx;i++){
      for(j=1;j<=Ly;j++){
	s = j+(i-1)*Lx+(k-1)*Lxy;
	if(particle[s] == 1){
	  fprintf(fp1,"1  1  1");
	  fprintf(fp1,"   ");
	};
	if(particle[s] == 2){
	  fprintf(fp1,"15  15  15");
	  fprintf(fp1,"   ");
	};
	if(particle[s] == 0){
	  fprintf(fp1,"9  9  9");
	  fprintf(fp1,"   ");
	};
      };
      fprintf(fp1,"\n");
    }; 
  };   
  fclose(fp1);      
}

///Subroutine to visualize the system///
void Visualization_povray(int tt, int ensemble){
  char fn[100];
  FILE *fp;
  
  int i,j,k,s,ncnt;
  int x,y,z;
  
  //printf("povray -D +A +J0 +Q10 %dsnp%dcnt_%dNch.pov\n",100000+tt,N_cent,N_ch);
  sprintf(fn,"%dsnp%dcnt_%dNch%dEnsemble.pov",100000+tt,N_cent, N_ch, ensemble);
  fp = fopen(fn,"w");
  
  fprintf(fp,"background{color rgb <1, 1, 1>}\n");
  fprintf(fp,"\n");
  
   /*********** Cell ***********/
  fprintf(fp,"sphere{<%d, %d, %d>, 1 scale<%d, %d, %d>\n",0,0,0,Ax,Ay,Az);
  //fprintf(fp,"translate<%d, %d, %d>\n",X0,Y0,Z0);
  fprintf(fp,"finish {phong 0.3}\n");
  fprintf(fp,"pigment{color rgbt<0.85, 0.7, 0.9, 0.8 >}}\n");
  fprintf(fp,"\n");
    /**************************/
  for(k=1;k<=Lz;k++){
    for(i=1;i<=Lx;i++){
      for(j=1;j<=Ly;j++){
	s = j + (i-1)*Lx + (k-1)*Lxy;
	
	//Translate the center to 0,0,0
	x = i-X0; y = j-Y0; z = k-Z0; 
	if(particle[s] == 1){
	  fprintf(fp,"sphere{<%d, %d, %d>, %f\n",x,y,z,0.5);
	  fprintf(fp,"pigment{color rgbt <0.9, 0.5, 0.5, 0.0>}\n");
	  fprintf(fp,"finish{ambient 0.5 specular 0.5}}\n");
	};
	if(particle[s] == 2){
	  fprintf(fp,"sphere{<%d, %d, %d>, %f\n",x,y,z,0.5);
	  fprintf(fp,"pigment{color rgbt <0.5, 0.9, 0.5, 0.0>}\n");
	  fprintf(fp,"finish{ambient 0.5 specular 0.5}}\n");
	};
      };
      //fprintf(fp,"\n");
    };    
  };
  
  
  
  //light source    
  fprintf(fp,"\n");
  fprintf(fp,"light_source\n");
  fprintf(fp,"{\n");
  
  fprintf(fp,"<%d, %d, %d> color rgb <1, 1, 1>\n",-(Lx+30),0,(Lz+30));
  fprintf(fp,"area_light <-150, 0, 0>, <0, 150, 0>, 25, 25\n");
  fprintf(fp,"jitter\n"); 
  fprintf(fp,"adaptive 1\n");
  fprintf(fp,"rotate <60, 60, 60>\n");
  fprintf(fp,"}\n");
  
  // camera position
  fprintf(fp,"\n");
  fprintf(fp,"camera\n");
  fprintf(fp,"{\n");
  fprintf(fp,"location <%d, %d, %d>\n",0,0,Lz+30);
  fprintf(fp,"look_at <%d, %d, %d>\n",0, 0,0);
  fprintf(fp,"rotate <0,90,0>");
  
  fprintf(fp,"\n");
  fprintf(fp,"}\n");
  fclose(fp);
  
}


///final positions of centrosomes and chromosomes////
void Final_configuration(int ensemble,int MAX_POWER){
  FILE *fp1, *fp2;
  char fn1[200],fn2[200];
  int i,j,k,s,x,y,z;
  
  sprintf(fn1,"Centrosome_final_config_time_power%d_Ncent%d_Nch%d_Fcs_csk%1.1f_Fcs_csd%1.1f_Fcs_kt%1.1f_Fcs_ch%1.1f_Fcs_surf%1.1f_ensemble%d_confined.dat",MAX_POWER,N_cent,N_ch,fcs_csk,fcs_csd,fcs_kt,fcs_ch,fcs_surf,ensemble);
  fp1=fopen(fn1,"w"); 
  sprintf(fn2,"Chromosome_final_config_time_power%d_Ncent%d_Nch%d_Fcs_csk%1.1f_Fcs_csd%1.1f_Fcs_kt%1.1f_Fcs_ch%1.1f_Fcs_surf%1.1f_ensemble%d_confined.dat",MAX_POWER,N_cent,N_ch,fcs_csk,fcs_csd,fcs_kt,fcs_ch,fcs_surf,ensemble);
  fp2=fopen(fn2,"w"); 
  
  for(k=1;k<=Lz;k++){
    for(i=1;i<=Lx;i++){
      for(j=1;j<=Ly;j++){
	s = j + (i-1)*Lx + (k-1)*Lxy;
        //Translate the center to 0,0,0
	x = i-X0; y = j-Y0; z = k-Z0; 
	if(particle[s]==1){
	  fprintf(fp1,"%d %d %d\n",x,y,z);
	}
	if(particle[s]==2){
	  fprintf(fp2,"%d %d %d\n",x,y,z);
	}
	
      }
    }
  }
  fclose(fp1);
  fclose(fp2);
}


/*************************/
void IdentifyPoles(){
  int i,j,k,l,p,q,r,sum,site,site1,site2;
  double dist;
  
  for(site=1;site<=N;site++){
    label_particle[site]=0;
    particle_dummy[site]=0;
  };
  
  for(k=1;k<=Lz;k++){
    for(i=1;i<=Lx;i++){
      for(j=1;j<=Ly;j++){
	site1 = j + (i-1)*Lx + (k-1)*Lxy;
	if((particle[site1]==1) && (label_particle[site1]==0)){//particle[site]=1 is centrosome
	  label_particle[site1]=1;
	  particle_dummy[site1]=1;//put a centrosome in the dummy lattice 
	  for(p=1;p<=Lz;p++){
	    for(q=1;q<=Lx;q++){
	      for(r=1;r<=Ly;r++){
		site2 = r + (q-1)*Lx + (p-1)*Lxy;
		if((particle[site2]==1) && (site2!=site1)){// && (label_particle[site2]==0)){
		  dist = sqrt((double)((p-k)*(p-k)+(q-i)*(q-i)+(r-j)*(r-j)));
		  if(dist<=dist_merge){//if centrosomes are closer than dist_merge, then merge them into one cent
		    if(label_particle[site2]==0)
		      label_particle[site2]=1; //include this particle into the same cluster if distance is less
		    else particle_dummy[site1]=0;
		  };
		};
	      };
	    };
	  };
	};
      };
    };
  };
  
  Npole_inst = 0;
  for(site=1;site<=N;site++){
    if(particle_dummy[site]==1){
      Npole_inst += 1;
      Npoles += 1;
    };
  };
  
}

////Subroutine to obtain the histogram of final number of centrosomal clusters
void OutputPolesHist(int ensemble){
  char fn[100];
  FILE *fp;
  
  int i,j,k,s,ncnt;
  int x,y,z;

  sprintf(fn,"HistPole%dcnt%dNch%1.2fMergeDist.dat",N_cent,N_ch,dist_merge);
  if(ensemble==1){
    fp = fopen(fn,"w");
  }else fp = fopen(fn,"a");
  
  fprintf(fp,"%d\n",Npole_inst);
  fclose(fp);
}


///Ensemble average of the final number of centrosomal clusters
void OutputPoles(){
  char fn[100];
  FILE *fp;
  
  int i,j,k,s,ncnt;
  int x,y,z;
  
  sprintf(fn,"Pole%dcnt%dNch%1.2fMergeDist.dat",N_cent,N_ch,dist_merge);
  fp = fopen(fn,"w");
  
  fprintf(fp,"%d %f\n",N_cent,(double)Npoles/(double)MAX_ENSEMBLE);
}


/***********************************************/
/*                                             */
/*             Mitotic Spindle                 */
/*      Centerosome-Chromosome Patterning      */
/*                                             */
/***********************************************/


main( int argc, char** argv )
{
  char fn[20], fn2[40], fn3[40], fn4[40];
  FILE *fp, *fp2, *fp3, *fp4;
  
  int i, j, k,r, x,y,n, idum, seed, ensemble, mcstp,mc_max;
  double p;
  
  /************************** PARAMETER TABLE ****************************/
  Lx = 60;  				    /// Linear system size along X-axis (variable)
  Ly = Lx;
  Lz = Lx;
  Lxy = Lx*Ly;
  MAX_POWER   = 13;                         /// Run time 2^13 (variable)
  
  N = Lx*Ly*Lz;                             /// Total number of lattice points
  N_cent = 8;                               ///Number of centrosomes
  N_ch = 46;                                ///Number of chromosomes
  mc_max = Pw(MAX_POWER);
  MIN_ENSEMBLE = 1;
  MAX_ENSEMBLE=200;                         ///No of ensemble
    
  N_total = N_cent+N_ch;                    ///Total number of particles
  Temp = 0.01;
  fcs_csd= -1.5;                            ///amplitude of inter-centrosomal attraction                       
  fcs_csk=  0.0;                     	    ///amplitude of inter-centrosomal repulsion   
  fcs_kt=-4.0;                              ///amplitude of centrosome-kinetochore attraction
  fcs_ch=10.0;                              ///amplitude of centrosome-chromosome repulsion
  fch_ch=4.0;                               ///steric repulsion between pairs of chromosomes
  fkt_kt=-0.0;                              ///long distance attraction between KTs
  r_cut = 2.0;                              ///cut-off distance for steric repulsion between pairs of chromosomes
  fch_wal = 2.0;                            ///steric repulsion between wall and chromosome
  fcs_wal= 0.0;                             ///steric repulsion between wall and centrosome
  fcs_surf = -0.4;                          ///amplitude of centrosome-cortex attraction
  
  //Center of the cell
  X0 = Lx/2; Y0 = Ly/2; Z0 = Lz/2;          ///center of the spheroidal cell
  Ax = Lx/3; Ay=Ly/4; Az=Lz/4;              ///semi-axes of the spheroidal cell 
  Ax_ini= Lx/4; Ay_ini= Ly/4; Az_ini= Lz/4; ///Defining initial volume for the distribution of centrosomes and chromosomes
  //Ax_ini=15; Ay_ini=15; Az_ini=15;
  Lav = 20.0;                               ///range of CS-CH and CS-CS interactions
  Lav_ast = Lav/4.0;                        ///range of CS-Cortex interaction
  
  dist_merge = 1.5;                         ///merging distance below which centrosomes are considered clustered
  /***********************************************************************/

  /********* Memory allocation ***********/
  occ                    = ivector (0,N); 
  particle   		 = ivector (0,N);
  neighb     		 = imatrix (0,N, 0,6); 
  label_particle 	 = ivector (0,N);
  particle_dummy 	 = ivector (0,N);
  
  
  Npoles=0;
  NearestNeighbAllocation();                                  ///Define nearest neighbor allocation
  
  for(ensemble = 1; ensemble<=MAX_ENSEMBLE; ensemble ++){     ///starting of ensemble loop 
    OccConfig();
    n =0 ;
    p = 0.0;
    //PrintPNGSnpshot(n);
    //t1 = Pw(n);
    t1 = Pw(p);
    //t1 = 5;
    
    for(mcstp=1; mcstp<=mc_max; mcstp++){                     ///starting of monte carlo step loop
      Monte_Carlo_Update();
      /* if( mcstp == t1){	
      //PrintPNGSnpshot(n);
      if (ensemble <= 20){
      //Visualization_povray(t1, ensemble);}
      //t1 +=50;
      if(mcstp<32){
      p += 1.0;
      }else{
      p += 0.1;
      };
      t1 = Pw(p);
      };*/
      //if(mcstp>=16000)fcs_cs=-4.0;
    };
    IdentifyPoles();
    OutputPolesHist(ensemble);
    Final_configuration(ensemble,MAX_POWER);
  }; 
  OutputPoles();
}
