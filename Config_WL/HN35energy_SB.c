/* gcc -O3 -lm  file.c -o  file.x  */
/* gcc -O3 -lm -lgsl -lgslcblas  file.c -o  file.x  */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fit.h>

/* #include<sys/types.h>
   #include<unistd.h>
   #include<term.h>
   #include<ncurses.h> */
#define min(A,B) ((A) < (B)) ? (A) : (B)
#define max(A,B) ((A) > (B)) ? (A) : (B)
#define MAXLINK 50     /* maximal number of links for a vertex */
#define BONDSCALE 1

typedef double TotalEnergy;
typedef double BondType;

struct Vertex {
  int links;               /* # of connections */
  BondType **bondcut; /* pointer to cutlist, which is +|bond| if bond 
				 is satisfied, -|bond| if bond is cut */
  struct Vertex **neighbor;   /* nearest neighbor connectivity vector of i-th edge to j neighbors */
  int index;   /* vertex' own index label, makes only sense for array */
};


gsl_rng *rnd;
int *minus;
int *vlist;
struct Vertex *V;
BondType *cutlist;
BondType *bondlist;
TotalEnergy energy_offset; /* offset only necessary for graph-redux */

BondType Babs(BondType x){
  if(x > 0){
    return(x);
  }else{
    return(-x);
  }
}
//argc,argv,&maxk,&k0,&k1,&y)
void commandlineparse(int argc, char **argv,int *k,BondType *k0,BondType *k1,BondType *y, int *seed){  
  int i;
  char outfile[100],mkdir[100];  
  FILE *fp;

  *seed = getpid();//default
  for (i = 1; i < argc; i++){  //Start at i = 1 to skip the command name.
    if (argv[i][0] == '-'){      switch (argv[i][1]){
      case 'k':       *k = atoi(argv[++i]);
        break;
/*       case 'i':       *maxinst = atoi(argv[++i]); */
/* 	break; */
/*       case 's':       *seed = atoi(argv[++i]); */
/*         break; */
      case '0':       *k0 = 0.01*atoi(argv[++i]);
        break;
      case '1':       *k1 = 0.01*atoi(argv[++i]);
        break;
      case 'y':       *y = 0.001*atoi(argv[++i]);
        break;
      default:
        fprintf(stderr,"\nError:  Incorrect option %s\n",argv[i]);
        fprintf(stderr,"\nAvailable options: \n\
  -k MaxSize of Hanoi Network N=2^k+1 \n\
  -0 100*K0, strength of BackBone Bonds \n\
  -1 100*K1, strength of SmallWorld Bonds \n\
  -y 1000*y, relative strength of L-Bonds to BackBone Bonds, L=y*K0 \n\
\n");
        exit(0);
      }//switch
    }//if
  }//for
  
/*   //set up output directory and log-file   */
/*   inttostring(*k,0,outfile); */
/*   fp = fopen(strcat(outfile,"/README"), "a"); */
/*   if(fp == NULL){ */
/*     inttostring(*k,0,outfile); */
/*     mkdir[0] = '\0'; */
/*     strcat(mkdir,"mkdir  "); */
/*     system(strcat(mkdir,outfile)); */
/*     fp = fopen(strcat(outfile,"/README"), "w"); */
/*     fprintf(fp,"N\tmaxinstance\tseed\tf\n"); */
/*   } */
/*   fprintf(fp,"%d\t%d\t%d\t%f\n",*k,*maxinst,*seed,*f); */
/*   fclose(fp); */
}




/******************************************************************/
/*                                                                */
/*      Input and output                                          */
/*                                                                */
/******************************************************************/

void opennewfile(char *fname,BondType k0,BondType k1,BondType y){
  int i=0;
  FILE *fp;

/*   i = inttostring(k,i,fname); */
/*   fname[i++] = '/'; */
  fname[i++] = 'H';
  fname[i++] = 'N';
  fname[i++] = '6';
  if(k0 < 0){
    fname[i++] = '-';
  }else{
    fname[i++] = '+';
  }
  if(k1 < 0){
    fname[i++] = '-';
  }else{
    fname[i++] = '+';
  }
  i = inttostring((int)(1000*y+.5),i,fname);
/*   fname[i++] = '.'; */
/*   i = inttostring(seed1,i,fname); */
  
  fp = fopen(fname,"w");
  fprintf(fp,"#\tk0\tk1\ty\n");
  fprintf(fp,"#\t%5.2f\t%5.2f\t%5.2f\t\n",(float)k0,(float)k1,(float)y);
  fclose(fp);
}


void out_run(char *fname,int k,TotalEnergy minenergy){
  FILE *fp;
  
  fp = fopen(fname,"a");
  fprintf(fp,"%4d  %20.10g\n",k,(double)minenergy/BONDSCALE);
  fclose(fp);
}

/*
  inputs an integer n and converts it into a sequence of characters in *str,
  starting at str[i] and returning the integer i' where "n" ends in *str 
  (str[i']='\0').
*/
int inttostring(int n,int i,char *str){
  if(n < 0){
    str[i++] = '-';
    i = inttostring(-n,i,str);
  }else{
    if((n/10) > 0) i = inttostring(n/10,i,str);
    str[i++] = '0' + n%10;
    str[i] = '\0';
  }
  return(i);
}

/********** END of Input/Output       *****************************/



/******************************************************************/
/*                                                                */
/*      Making of a random graphs                                 */
/*                                                                */
/******************************************************************/


/*
int initHanoiGraphN5(int nv,float perc){
  int numlines=0,l;
  int i,j,k,dx;
  
  for(i = 0; i<=nv; i++){
    V[i].index = i;
    V[i].links = 0;
  }
  //HN3 links
  k = nv;
  while( (k/=2) > 1 ){
    i = k/2;
    while( i<nv ){
      j = i+k;
      //   printf(" %d %d\n",i,j);
      if(gsl_rng_uniform(rnd) < perc){
	V[i].neighbor[V[i].links++] = &V[j];
	V[j].neighbor[V[j].links++] = &V[i];
	++numlines;
      }
      i = j+k;
    }
  }

  // HN5 links 
  k = 2*nv;
  while( (k/=2) > 1 ){
    i = 0;
    while( i<nv ){
      j = i+k;
      //     printf(" %d %d\n",i,j);
      if(gsl_rng_uniform(rnd) < perc){
	V[i].neighbor[V[i].links++] = &V[j];
	V[j].neighbor[V[j].links++] = &V[i];
	++numlines;
      }
      i = j;
    }
  }
  for(i=0; i<nv; i++){
    j = i+1;
    if(j <= nv){
      //     printf(" %d %d\n",i,j);
      if(gsl_rng_uniform(rnd) < perc){
	V[i].neighbor[V[i].links++] = &V[j];
	V[j].neighbor[V[j].links++] = &V[i];
	++numlines;
      }
    }
  }
  // printf("\n");
  return(numlines);
}
*/

//initHanoi5graph(nv,k0,k1,y)
int initHanoi5graph(int nv,BondType k0, BondType k1,BondType y){
  int l,k,flag;
  int ii,i,j;
  
  for(i = 0; i<=nv; i++){
    V[i].index = i;
    V[i].links = 0;
  }
  l = 0;
  k = 2;
  do{
    flag = 1;
    for(i = k/2 ; i<nv; i+=k){
      if(flag){//make HN3 bonds
	j = i+k;
	//        printf(" %5d %5d %5d\n",l,i,j);
	cutlist[l] = k1*BONDSCALE; //small-world bond        
	V[i].bondcut[V[i].links] = &cutlist[l];
	V[j].bondcut[V[j].links] = &cutlist[l];
	V[i].neighbor[V[i].links++] = &V[j];
	V[j].neighbor[V[j].links++] = &V[i];
	++l;
      }
      flag = !flag;

//this adds HN5-bonds, get HN3 for y=0!
      ii= i-k/2;
      j = i+k/2;
      //        printf(" %5d %5d %5d\n",l,i,j);
      cutlist[l] = y*k0*BONDSCALE; //L-bond      
      V[ii].bondcut[V[ii].links] = &cutlist[l];
      V[j].bondcut[V[j].links] = &cutlist[l];
      V[ii].neighbor[V[ii].links++] = &V[j];
      V[j].neighbor[V[j].links++] = &V[ii];
      ++l;
    }
  }while( (k *= 2) < nv);
  
  for(i = 0; i<nv; i++){
    j = i+1;
    //  printf(" %5d %5d %5d\n",l,i,j);
    cutlist[l] = k0*BONDSCALE; //backbone bond        
    V[i].bondcut[V[i].links] = &cutlist[l];
    V[j].bondcut[V[j].links] = &cutlist[l];
    V[i].neighbor[V[i].links++] = &V[j];
    V[j].neighbor[V[j].links++] = &V[i];
    ++l;
  }
  return(l);
}




int initHanoi6graph(int nv,BondType k0, BondType k1,BondType y){
  int l,k,flag;
  int ii,i,j;
  
  for(i = 0; i<=nv; i++){
    V[i].index = i;
    V[i].links = 0;
  }
  l = 0;
  k = 2;
  do{
    flag = 1;
    for(i = k/2 ; i<nv; i+=k){
      j = i+flag*(k+k/2);
      //        printf(" %5d %5d %5d\n",l,i,j);
      cutlist[l] = k1*BONDSCALE; //small-world bond        
      V[i].bondcut[V[i].links] = &cutlist[l];
      V[j].bondcut[V[j].links] = &cutlist[l];
      V[i].neighbor[V[i].links++] = &V[j];
      V[j].neighbor[V[j].links++] = &V[i];
      flag = -flag;
      ++l;
      //this adds HN6-bonds, get HNNP for y=0!
      ii= i-k/2;
      j = i+k/2;
      //        printf(" %5d %5d %5d\n",l,i,j);
      cutlist[l] = y*k0*BONDSCALE; //L-bond      
      V[ii].bondcut[V[ii].links] = &cutlist[l];
      V[j].bondcut[V[j].links] = &cutlist[l];
      V[ii].neighbor[V[ii].links++] = &V[j];
      V[j].neighbor[V[j].links++] = &V[ii];
      ++l;
    }
  }while( (k *= 2) < nv);
  
  for(i = 0; i<nv; i++){
    j = i+1;
    //  printf(" %5d %5d %5d\n",l,i,j);
    cutlist[l] = k0*BONDSCALE; //backbone bond        
    V[i].bondcut[V[i].links] = &cutlist[l];
    V[j].bondcut[V[j].links] = &cutlist[l];
    V[i].neighbor[V[i].links++] = &V[j];
    V[j].neighbor[V[j].links++] = &V[i];
    ++l;
  }
  return(l);
}



/******************************************************************/
/*                                                                */
/*      Fixing bonds and initializing spins                       */
/*                                                                */
/* "initbonds" makes "bondlist" with +/-1, "initspins" assigns    */
/* +/-1 spins to each vertex and marks bonds either violated (cut)*/
/* or not (state -1 or +1) in array "cutlist".                    */
/* IMPORTANT: "V[n].bondcut[i]" points to respective element in   */
/* in "cutlist" (not "bondlist"!), ie. it does NOT refer to the   */
/* weight of the bond but its STATE! Thus, in an update (see      */
/* "flip") it is easier to determine all fitnesses.               */
/*                                                                */
/******************************************************************/

/*
TotalEnergy initGauss(int numlines){
  int i;
  TotalEnergy balance=0;
  
  for(i=0; i<numlines; i++){
    bondlist[i] = cutlist[i] =  (BondType)gsl_ran_gaussian(rnd,BONDSCALE);
    balance += Babs(cutlist[i]);
  }
  return(balance);
}
*/

/*
  "fixdoublebond" checks whether two vertices v1 and v2 are doubly connected.
  If their are, their double-link can be replaced by a single link according to
  new_j12=j12.1+j12.2, which may be zero in which case j1 and j2 become 
  disconnected. In either case, we have to consider the effect for the overall 
  cost: if both bonds are equal their cost is either zero or their sum, if both
  bonds are opposite, one of them is always cut!
*/



void removelink(struct Vertex *v1, struct Vertex *v2, int i){
  int j=0;    
  
  *(v1->bondcut[i]) = 0;   /* forget the second bond */
  while(v2->bondcut[j] != v1->bondcut[i]) ++j;
  v2->bondcut[j] = v2->bondcut[--(v2->links)];        /* override link... */
  v2->neighbor[j] = v2->neighbor[v2->links];    /* ...with new neighbor */
  v1->bondcut[i] = v1->bondcut[--(v1->links)];     /* override link... */
  v1->neighbor[i] = v1->neighbor[v1->links]; /* ...with new neighbor */
}


void fixdoublebond(struct Vertex *v1, struct Vertex *v2){
  int i[20],j;
  int k;
  
 restart:
  k=0;
  for(j=(v1->links)-1; j>=0; j--){
    if(v1->neighbor[j] == v2) i[k++] = j;
  }
  if(k >= 2){
    //    printf("error in fixdoublebond: %d %d\n",v1->index,v2->index); /* at most have a double bond! */
    //  }else if(k == 2){
    *(v1->bondcut[i[0]]) += *(v1->bondcut[i[1]]);/* j12.1=j12.1+j12.2 */
    if(*(v1->bondcut[i[0]]) == 0) removelink(v1,v2,i[0]);
    removelink(v1,v2,i[1]);
  }else if(k == 1){
    if(*(v1->bondcut[i[0]]) == 0) removelink(v1,v2,i[0]);
  }
  if(k >2) goto restart;
}

/* /\* */
/* "graph_redux" eliminates all 0-,1-,2-,and 3-connected vertices from a graph. */
/* *\/ */

/* void graph_redux(int *nv,int numlines){   */
/*   int len,l,smallest_link,current,previous,n,i,j; */
/*   //     testgraph(*nv); */
/*   for(i=0; i<*nv; i++) */
/*     for(j=0; j< *nv; j++)fixdoublebond(&V[i],&V[j]); */
/*        testgraph(*nv); */
/*   for(l = 0; l<numlines; l++) cutlist[l] = bondlist[l]; /\*sync cutlist *\/ */
/*   smallest_link = MAXLINK; */
/*   current = 0; */
/*   len = 0; */
/*   for(n=0; n<*nv; n++) */
/*     if(V[n].links > 0){  //ignore 0-connected vertices */
/*       vlist[current] = n; */
/*       current = n; */
/*       ++len; */
/*       smallest_link = min(smallest_link,V[n].links); */
/*     } */
/*   vlist[current] = ((V[0].links > 0)?0:vlist[0]);//close the loop */
  
/*   while(smallest_link < 4){ */
/*   restart: */
/*     for(l=0; l<len; l++){ */
/*       current = vlist[previous=current]; */
/*       //      printf("\n%d %d >",smallest_link,n=current); */
/*       //      for(i=0; i<len; i++)printf("%4d",n=vlist[n]); */
/*       //      printf("\n"); */
/*       if((V[current].links <= smallest_link)){// || superbond(&V[current])){ */
/* 	n = reduce_vertex(&V[current]); //eliminate current vertex */
/* 	if(testgraph(*nv)) printf("%d: %d %d %d %20.10g\n",n,smallest_link,current,len,energy_offset);  */
/* 	vlist[previous] = vlist[current]; //forget current vertex */
/* 	current = previous; //start over counting */
/* 	smallest_link = 1; //may have new 0- and 1-connects */
/* 	--len; */
/* 	goto restart; */
/*       }//end if V<=smallest */
/*     }//end for */
/*     ++smallest_link; */
/*   } */
  
/*   for(l = 0; l<numlines; l++) bondlist[l] = cutlist[l];/\*re-sync bondlist *\/ */
/*   //  reform_graph(nv); */
/*   //  printf("=========================================\n"); */
/* } */

int reduce_vertex(struct Vertex *v){
  //     printf("==================== %1d-redux on %3d-vertex\n",v->links,v->index);
  switch(v->links){
  case 0: //ignore and forget
    return(0);
    break;
  case 1: one_redux(v);
    //    printf("1  ");
    return(1);
    break;
  case 2: two_redux(v);
    //   printf("2   ");
    return(2);
    break;
  case 3: three_redux(v);
    //  printf("3   ");
    return(3);
    break;
  default: 
    return(v->links);
  }//end switch
}

/* void reform_graph(int *nv){ */
/*   int n,i,j; */
/*   struct Vertex *v; */
  
/*   n = 0; */
/*   while((*nv > 0) && (V[*nv-1].links <= 0)) --(*nv); */
/*   while(n < *nv){ */
/*     if(V[n].links <= 0){ */
/*       --(*nv); */
/*       V[n].links = V[*nv].links; */
/*       for(i=0; i<V[*nv].links; i++){ */
/* 	V[n].bondcut[i] = V[*nv].bondcut[i]; */
/* 	V[n].neighbor[i] = V[*nv].neighbor[i]; */
/* 	v = V[*nv].neighbor[i]; */
/* 	for(j=0; j<(v->links); j++){ */
/* 	  if(v->neighbor[j] == &V[*nv]) v->neighbor[j] = &V[n]; */
/* 	} */
/*       } */
/*       while((*nv > 0) && (V[*nv-1].links <= 0)) --(*nv); */
/*     } */
/*     ++n; */
/*   } */
/* } */



/*
  "one_redux" sets the #-of-links of v to zero and unlinks the edge to v in
  the neighbor w of v. Note that the pointer to the respective bond gets 
  overridden, but there is no need to eliminate it from "bondlist".
*/ 
int one_redux(struct Vertex *v){/* vertex v has only one neighbor! */
  struct Vertex *w;
  int i=0;
  
  w = v->neighbor[0];    /* fix neighboring vertex w */
  while(w->neighbor[i] != v) ++i; 
  w->neighbor[i] = w->neighbor[--(w->links)];
  w->bondcut[i] = w->bondcut[w->links];
  v->links = -1;  /*  eliminate vertex v; need -1 instead of 0; note that it is important to allow w->links to be zero. Only then the degeneracy is calculated correctly in "reform_graph".  */
  energy_offset -= (TotalEnergy)Babs(*(v->bondcut[0])); /* E = E' - |bond_to_v| */
  *(v->bondcut[0]) = 0;     /* set respective bondweight zero */
  return(w->index);
}



/*   
     "two_redux" eliminates a two-connected vertex v by reconnecting the adjacent 
     vertices w1 and w2 with the appropriate reduced bond. The new bond gets the 
     value j12=(|j1+j2|-|j1-j2|)/2 [and the energy is offset by 
     de=(|j1+j2|-|j1-j2|)/2], since for optimally chosen v, 
     -(j1*v*w1+j2*v*w2)=-|j1*w1+j2*w2|=-(de+j12*w1*w2)
     where j1 linked v to w1 and j2 linked v to w2.
     Important: we also have to check whether w1 and w2 now have a double bond!
*/

int two_redux(struct Vertex *v){/* v has exactly two neighbors! */
  struct Vertex *w1,*w2;
  int i1=0,i2=0;
  BondType J12,J1,J2;
  
  w1 = v->neighbor[0];
  w2 = v->neighbor[1];
  while(w1->neighbor[i1] != v) ++i1; 
  w1->neighbor[i1] = w2;
  while(w2->neighbor[i2] != v) ++i2; 
  w2->neighbor[i2] = w1;
  
  J1 = *(w1->bondcut[i1]);   /*  Reconnect w1, w2 with new bond J12 */
  J2 = *(w2->bondcut[i2]);
  J12 = (Babs(J1+J2)-Babs(J1-J2))/2;
  *(w1->bondcut[i1]) = J12;
  *(w2->bondcut[i2]) = 0;    /*  Forget old bond  */
  w2->bondcut[i2] = w1->bondcut[i1];  /* redirect bond */
  fixdoublebond(w1,w2);      /* fix if w1,w2 are doubly connected now */
  v->links = -1;                  /*  eliminate vertex v; need -1 instead of 0; note that it is important to allow w->links to be zero. Only then the degeneracy is calculated correctly in "reform_graph".  */
  energy_offset -= (TotalEnergy)(Babs(J1+J2)+Babs(J1-J2))/2; 
  return(min(w1->index,w2->index));
}



int three_redux(struct Vertex *v){
  struct Vertex *w1,*w2,*w3;
  int i1=0,i2=0,i3=0;
  BondType J12,J13,J23,w,x,y,z,J1,J2,J3;
  
  J1 = *(v->bondcut[0]);  
  J2 = *(v->bondcut[1]);
  J3 = *(v->bondcut[2]);
  w = Babs(J1-J2+J3); 
  x = Babs(-J1+J2+J3); 
  y = Babs(J1+J2+J3); 
  z = Babs(-J1-J2+J3);
  J12 = (-w-x+y+z)/4;
  J23 = (-w+x+y-z)/4;
  J13 = (w-x+y-z)/4;
  *(v->bondcut[0]) = J23;  
  *(v->bondcut[1]) = J13;
  *(v->bondcut[2]) = J12;
  
  w1 = v->neighbor[0];
  w2 = v->neighbor[1];
  w3 = v->neighbor[2];
  while(w1->neighbor[i1] != v) ++i1; 
  w1->neighbor[i1] = w2;
  w1->neighbor[w1->links] = w3;
  while(w2->neighbor[i2] != v) ++i2; 
  w2->neighbor[i2] = w3;
  w2->neighbor[w2->links] = w1;
  while(w3->neighbor[i3] != v) ++i3; 
  w3->neighbor[i3] = w1;
  w3->neighbor[w3->links] = w2;
  
  w1->bondcut[i1] = v->bondcut[2];         /* Let new w1-w2 bond use...*/
  w2->bondcut[w2->links++] = v->bondcut[2];/* ...old v-w3 bond */
  w2->bondcut[i2] = v->bondcut[0];         /* Let new w2-w3 bond use...*/
  w3->bondcut[w3->links++] = v->bondcut[0];/* ...old v-w1 bond */
  w3->bondcut[i3] = v->bondcut[1];         /* Let new w1-w3 bond use...*/
  w1->bondcut[w1->links++] = v->bondcut[1];/* ...old v-w2 bond */
  
  fixdoublebond(w1,w2);      /* fix if w1,w2 are doubly connected now */
  fixdoublebond(w1,w3);  
  fixdoublebond(w2,w3);  
  v->links = 0;                  /*  eliminate vertex v  */
  energy_offset -= (TotalEnergy)(w+x+y+z)/4;
  return(min(min(w1->index,w2->index),w3->index));
}




/* int testgraph(int nv){ */
/*   struct Vertex *v; */
/*   int flag,i,j,k,error=1; */
  
/*   for(i=0; i<nv; i++){ */
/*     for(j=0; j<V[i].links; j++){ */
/*       v=V[i].neighbor[j]; */
/*       flag=0; */
/*       for(k=0; k<(v->links); k++) if(v->neighbor[k] == &V[i]) flag=1; */
/*       if(!flag){ */
/* 	printf("Trouble with %d %d\n",i,v->index); */
/* 	error=1; */
/*       } */
/*     } */
/*   } */
/*   if(error){ */
/*     printf("\n"); */
/*     for(i=0; i<nv; i++) if(V[i].links>0) printf("%4d",i); */
/*     printf("\n"); */
/*     //    for(i=0; i<nv; i++) if(V[i].links>0) printf("%4d",V[i].slaves); */
/*     //    printf("\n"); */
/*     for(i=0; i<nv; i++) if(V[i].links>0) printf("%4d",V[i].links); */
/*     printf("\n"); */
/*     flag=1; */
/*     j=0; */
/*     while(flag){ */
/*       flag=0; */
/*       for(i=0; i<nv; i++) { */
/* 	if(V[i].links>j){ */
/* 	  flag=1; */
/* 	  printf("%4d",(V[i].neighbor[j])->index); */
/* 	}else{ */
/* 	  if(V[i].links>0) printf("    "); */
/* 	} */
/*       } */
/*       printf("\n"); */
/*       ++j; */
/*     } */
/*     printf("\n"); */
/*     for(i=0; i<nv; i++) for(j=0; j<V[i].links; j++) if(i < (V[i].neighbor[j])->index) printf("%4d",i); */
/*     printf("\n"); */
/*     for(i=0; i<nv; i++) for(j=0; j<V[i].links; j++) if(i < (V[i].neighbor[j])->index) printf("%4d",(V[i].neighbor[j])->index); */
/*     printf("\n"); */
/*     for(k=0; k<1; k++){ */
/*       for(i=0; i<nv; i++) for(j=0; j<V[i].links; j++) if(i < (V[i].neighbor[j])->index) if(V[i].links>0) printf("%4.1f",(float)(*(V[i].bondcut[j]+k))); */
/*       printf("\n"); */
/*     } */
/*   } */
/*   return(error); */
/* } */


int testgraph(int nv){
  struct Vertex *v;
  int flag,i,j,k,error=1;
  
  for(i=0; i<nv; i++){
    for(j=0; j<V[i].links; j++){
      v=V[i].neighbor[j];
      flag=0;
      for(k=0; k<(v->links); k++) if(v->neighbor[k] == &V[i]) flag=1;
      if(!flag){
	printf("Trouble with %d %d\n",i,v->index);
	error=1;
      }
    }
  }
  if(error){
    printf("\n");
    for(i=0; i<nv; i++) if(V[i].links>0) printf("%4d",i);
    printf("\n");
    //    for(i=0; i<nv; i++) if(V[i].links>0) printf("%4d",V[i].slaves);
    //    printf("\n");
    for(i=0; i<nv; i++) if(V[i].links>0) printf("%4d",V[i].links);
    printf("\n");
    flag=1;
    j=0;
    while(flag){
      flag=0;
      for(i=0; i<nv; i++) {
	if(V[i].links>j){
	  flag=1;
	  printf("%4d",(V[i].neighbor[j])->index);
	}else{
	  if(V[i].links>0) printf("    ");
	}
      }
      printf("\n");
      ++j;
    }
    printf("\n");
    for(i=0; i<nv; i++) for(j=0; j<V[i].links; j++) if(i < (V[i].neighbor[j])->index) printf("%4d",i);
    printf("\n");
    for(i=0; i<nv; i++) for(j=0; j<V[i].links; j++) if(i < (V[i].neighbor[j])->index) printf("%4d",(V[i].neighbor[j])->index);
    printf("\n");
    for(k=0; k<1; k++){
      for(i=0; i<nv; i++) for(j=0; j<V[i].links; j++) if(i < (V[i].neighbor[j])->index) if(V[i].links>0) printf("%4.1f",(float)(*(V[i].bondcut[j]+k)));
      printf("\n");
    }
  }
  return(error);
}




/*************************************************************/
/*                                                           */
/*                                                           */
/*                 Main Program                              */
/*                                                           */
/*                                                           */
/*                                                           */
/*************************************************************/

int main(int argc, char **argv){ 
  int numlines,k,i,j,maxbonds;
  int seed,inst,repeat;
  TotalEnergy balance,defect_energy[2];
  int nv,nnv,maxk;
  BondType k0,k1,y;
  char outfile[100];

/*   rnd = gsl_rng_alloc(gsl_rng_mt19937); */
  commandlineparse(argc,argv,&maxk,&k0,&k1,&y, &seed);
  opennewfile(outfile,k0,k1,y);
    nnv = ldexp(1,maxk);  /* =2^k !? */     
    maxbonds = (int)(3.5*nnv);
    minus = (int *)malloc(maxbonds * sizeof(int));
    bondlist = (BondType *)malloc(maxbonds * sizeof(BondType));
    cutlist = (BondType *)malloc(maxbonds * sizeof(BondType));
    V = (struct Vertex *)malloc((nnv+1) * sizeof(struct Vertex));
    vlist = (int *)malloc((nnv+1) * sizeof(int));
    j = 1;
    while( (j *= 2) <= nnv){
      for(i = j/2 ; i<nnv; i+=j){
	V[i].bondcut = (BondType **)malloc(3*j * sizeof(BondType));
	V[i].neighbor = ((struct Vertex **)malloc(3*j * sizeof(struct Vertex)));
      }
    }
  V[0].bondcut = (BondType **)malloc(3*j * sizeof(BondType));
  V[0].neighbor = ((struct Vertex **)malloc(3*j * sizeof(struct Vertex)));
  V[nnv].bondcut = (BondType **)malloc(3*j * sizeof(BondType));
  V[nnv].neighbor = ((struct Vertex **)malloc(3*j * sizeof(struct Vertex)));

 
  for(k=2; k<maxk; k++){
    energy_offset = 0;
    nnv = ldexp(1,k);  /* =2^k !? */ 
    nv = nnv;
    /*     gsl_rng_set(rnd, ((unsigned long int) abs(((inst+1)*k*seed)))); */
    /*     for(j = 0; j <10; j++) { */
    /*       gsl_rng_uniform(rnd); //Do some initial cycling on the random generator */
    /*     } */
    numlines = initHanoi5graph(nv,k0,k1,y);//get HNNP for y=0
    j = 1;
    while( (j *= 2) <= nv){
      //            testgraph(nv+1);
      for(i = j/2 ; i<nv; i+=j) reduce_vertex(&V[i]);
    }
    //    testgraph(nv+1);
    reduce_vertex(&V[0]);
    out_run(outfile,k,energy_offset);
  }
  free(bondlist);
  free(cutlist);
  free(V);
  free(minus);
  free(vlist);
  return(0);
}

/********** END of MAIN           *** *****************************/
