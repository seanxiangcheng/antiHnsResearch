 /*
  gcc Program.c -lgsl -lgslcblas -lm -o program.x
	
  This code is to use Wang-Landau sampling to
  find the density of states of the lattice gas problem 
	with at-most-1-neighbor density constraint in HN3; 
	Periodic Boundary Condition is used in this code; 
	The output is a file with number of occupied particles and their 
 	corresponding density of states
  
  There are 2 other versions of this method for HN3 and HN5, 
  in which the data structure is a little different.
  
  For more information of the model and method:
  http://arxiv.org/abs/1409.8313
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <term.h>
#include <ncurses.h>
#include <gsl/gsl_rng.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>

#define MAXNEIG 36
#define MIN_INV_MU 0.05
#define INV_MU_STEP 0.001
#define INV_MU_L 800

int L=0, K=0, Flat=80, Seed, Name = 0; 	// Input parameters
float Flatness=0.8;
int StatesLength = 0, MaxLength, MaxState, CurrentState, CurrentIndex, MinState;
double Logfnow=1.0, Logfmin=0.00000001, MinStateLogg, MaxSteps=2.0e13;
double Initial_Steps=0.0, Check_Every_Steps=0.0, Total_Steps=0.0, Inner_Loop_Steps=0.0; 

gsl_rng *rnd;	// RNG: random number

struct StateInfo{
  int E;
  float H; 
  double logg;
}; 


/* 	integer power of an integer (with 0^0 = 1) 	*/
int Intpow(int base, int power){
  int x=1;
  int i;
  for(i = 0; i<power; i++) {x=x*base;}
  return(x);
}


/* 	integer log2 of an integer	*/
int Intlog2(int num){
  int power=0;
  while(num>1){
    num=num/2;
    power++;
  }
  return(power);
}


/* 	Find the level i of site n and sequential number j	*/
int *Findij(int n){
  static int ij[2];
  ij[0]=1;
  ij[1]=0;
  if(n==0)
  {
      ij[0] = 0;
      return(ij);
  }
  while(n%2==0) {
    n=n/2;
    ij[0]++;
  }
  ij[1]=(n-1)/2;
  return(ij);
}


/* Find the neighbors of site n in the complex network*/
int * FindNeighbors(n){ 
  	static int neighbors[MAXNEIG];
  	int *ij;
  	unsigned int i,index;
  	int ni0, nj0, ni1,ni2, n1,n2;
  	if(n>L){
    		fprintf(stderr,"\n Error:  site %d is out of the network length  %d\n", n,L);
    		return(0);
  	}
  	for(i=0; i<MAXNEIG; i++){
    		neighbors[i]=-1; 	// Fill up with real neighbors later so that we know where real neighbors ends
  	}
  	// initialize neighs first to be all 0
  	ij=Findij(n);
  	ni0=ij[0];
  	nj0=ij[1];

    if(n==0){
        neighbors[0]=1;
        index=1;
        for(i=0; i<(K-1); i++){
            neighbors[index++] = Intpow(2, i)*3;
        }
        return(neighbors);
    }
  	else if(ni0==1){	//  if level 1
    		if(nj0==0){
            neighbors[0]=0;
      			neighbors[1]=2;
      			neighbors[2]=4;
    		}
    		else if(nj0%2==0){
      			neighbors[0]=n-1;
      			neighbors[1]=n+1;
      			neighbors[2]=n+3;  
    		}
    		else{
      			neighbors[0]=n-3;
      			neighbors[1]=n-1;
      			neighbors[2]=n+1;
	    	}
    		return(neighbors);   
  	}
  	else if(ni0==2){	//  if level 1
    		if(nj0%2==0){
      			neighbors[0]=n-1;
      			neighbors[1]=n+1;
      			neighbors[2]=n+6;  
    		}
    		else{
      			neighbors[0]=n-6;
      			neighbors[1]=n-1;
      			neighbors[2]=n+1;
	    	}
    		return(neighbors);   
  	}
  	else if(ni0<K){
    		neighbors[0]=n-1;
    		neighbors[1]=n+1;
    		if(nj0%2==0){
      			neighbors[2]=n+3*Intpow(2,ni0-1);
    		}
    		else{
      			neighbors[2]=n-3*Intpow(2,ni0-1);
    		}
    		index=3;
    		if(ni0>2){
      			for(i=0;i<(ni0-2);i++){	// sites nn that may connect to site n by connecting back and forth; 
        			n1=n-3*Intpow(2,i);
        			n2=n+3*Intpow(2,i);
        
        			ij=Findij(n1);
        			ni1=ij[0];

        			ij=Findij(n2);
        			ni2=ij[0];

        			if(ni1==i+1){
          				neighbors[index++]=n1;
        			}
        
        			if(ni2==i+1){
          				neighbors[index++]=n2;
        			}
      			}  
    		}
    		return(neighbors);  
  	}
  	else if(ni0==K){
    		neighbors[0]=n-1;
    		neighbors[1]=n+1;
    		index=2;
    		for(i=0;i<(ni0-2);i++){
        		n1=n-3*Intpow(2,i);
        		n2=n+3*Intpow(2,i);
        		ij=Findij(n1);
        		ni1=ij[0];
        
        		ij=Findij(n2);
        		ni2=ij[0];

        		if(ni1==i+1){
          			neighbors[index++]=n1;
        		}
        
        		if(ni2==i+1){
          			neighbors[index++]=n2;
        		}
      		}      
  	}
  	else if(ni0==K+1){
    		neighbors[0]=L-1;
    		index=1;
    		for(i=0;i<(ni0-2);i++){
        		n1=n-3*Intpow(2,i);
        		ij=Findij(n1);
        		ni1=ij[0];
        		if(ni1==i+1){
          			neighbors[index++]=n1;
        		}
     		}          
  	}
  	else{
			fprintf(stderr,"\n Error in function FindNeighbors(): neighbors[3]<0. n=%d, L=%d\n", n,L);
  	} 
  	return(neighbors);
}



/*	Print the state.E, state.logg, and state.H for all the states found */
int Print_States(struct StateInfo state[]){
    int i;
    printf(" Print_States:\n");
    printf(" Current State: %d;\n", CurrentState);
    printf(" StatesLength: %d; \n", StatesLength);
    printf(" E          logg           H\n");
    for(i = 0; i < StatesLength; i++){
        printf(" %-12d %-18.10e %-12.2e\n", state[i].E, state[i].logg, state[i].H);
    }
    return(0);
}


/*	Initialize the StateInfo, neighbors and other parameters*/
void Initialize(struct StateInfo state[], int site[], int neighbors[], int neighloc[]){
  	int i, j;
  	int *neigh, neig_index = 0;
    for(i=0; i< (L+1); i++){
      site[i] = 1;
      neigh=FindNeighbors(i);
      neighloc[i]=neig_index;
      for(j=0; neigh[j]!=-1; j++){
        neighbors[neig_index++]=neigh[j];
      }
    }
    neighloc[L+1]=neig_index;
    printf("\n\n Total neig links: %d \n\n", neig_index);
    state[0].E = (Intpow(2,K+2)-4)/2;
    state[0].logg = 0.0;
    state[0].H = 0.0;
  	for(i = 1; i < MaxLength; i++)
    {
        state[i].E = -1;
        state[i].logg = 0.0;
        state[i].H = 0.0;
  	}
    StatesLength = 1;
    CurrentState = (Intpow(2,K+2)-4)/2;
    CurrentIndex = 0;
    MaxState = CurrentState;
    MinState = CurrentState;
}

/* Find the energy change if exchange size i and h */
int Exchange_Energy(int site[], int neighbors[], int i, int h){
  int j, e1 = 0, e2 = 0, e, neig_flag = 0;
  if(site[i] == site[h]){
    fprintf(stderr, "\n Error: In Exchange_Energy, site[%d] = %d, site[%d] = %d;\
     \n No need to exchange!", i, site[i], h, site[h]);
  }
  for(j = 0; j < 3; j++){
      if(neighbors[3*i+j] == h){
        neig_flag = 1;
        break;
      }
  }
  if(neig_flag){
    for(j = 0; j < 3; j++){
        e1 += site[neighbors[3*i + j]];
        e2 += site[neighbors[3*h + j]];
    }
    e = CurrentState + (site[h] - site[i])*(e1 - site[h] - e2 + site[i]);
    return(e);
  }
  else{
    for(j = 0; j < 3; j++){
      e1 += site[neighbors[3*i + j]];
      e2 += site[neighbors[3*h + j]];
    }
    e = CurrentState - 2*e1*site[i] - 2*e2*site[h];
  }
  return(e);
}

/* Binary Search based on sorted state.E */
int Binary_Search(struct StateInfo state[], int e, int imin, int imax){
  int imid;
  while(imin <= imax){
    imid = (imin + imax) / 2;
    if(state[imid].E == e) 
      return(imid);
    else if(state[imid].E < e) 
      imax = imid - 1;
    else 
      imin = imid + 1;
  }
  return(-1); // if nothing is found in the loop
}

/* Find the index in the most simple and cost way before sorting */
int Find_E_Index(struct StateInfo state[], int e){
  int index = 0;
  if(Total_Steps < Initial_Steps){
    for(index = 0; index < StatesLength; index++){
      if(state[index].E == e){
          return(index);
          break;
      }
    }
    if(index == StatesLength){//could remove this if, but keep the statements in it
      StatesLength++;
      index = StatesLength - 1;
      if(state[index].E != -1){ // could delete
        fprintf(stderr, "\n Error: The new state's state[%d].E != -1 \n\
        (all state.E is initialized to be -1)!!!\n StatesLength = %d\n \
             Current state[%d].E = %d", index, StatesLength-1, index, state[index].E);
        Print_States(state);
        exit(0);
      } 
      state[index].E = e; // logg, H is already initialized in function Initialize
      return(index);
    }
    else//could delete
      fprintf(stderr, "\n Error: In intial steps, StatesLength = %d, Index = %d", StatesLength, index); 
  }
  else{ // find sorted state in log(N) steps
    index = Binary_Search(state, e, 0, StatesLength - 1);
    if(index != -1)
      return(index);
    else
      return(-1);
  }
  return(index);
}

/*	Update the system by adding or removing particles*/
int Update(struct StateInfo state[], int site[], int neighbors[], int neighloc[]){
	int i, j, e2 = 0, ie1, ie2, begin, end;
  Total_Steps += 1.0;
	i = (int)((L+1)*gsl_rng_uniform(rnd));
  begin = neighloc[i];
  end = neighloc[i+1];
  for(j = begin; j < end; j++){
    e2 += site[neighbors[j]];
  }
  e2 = CurrentState - 2 * e2 * site[i];
  ie1 = CurrentIndex;
  ie2 = Find_E_Index(state, e2);
  if(ie2 == -1){
      site[i] = -site[i];
      CurrentState = e2;
      CurrentIndex = StatesLength;
      if(CurrentState < MinState)
        MinState = CurrentState;	
      return(e2);
  }
	if(state[ie1].logg > state[ie2].logg || exp(state[ie1].logg - state[ie2].logg) > gsl_rng_uniform(rnd)){
		site[i] = -site[i];
		CurrentState = e2;
    CurrentIndex = ie2;
		state[ie2].logg += Logfnow;
		state[ie2].H += 1.0;
    if(CurrentState < MinState)
      MinState = CurrentState;		
  }
  else{
    state[ie1].logg += Logfnow;
		state[ie1].H += 1.0;
	}
  /* Exchanging two sites */
  /*
	i = (int)(L * gsl_rng_uniform(rnd));
  h = (int)(L * gsl_rng_uniform(rnd));
  while(h == i){
    h = (int)(L * gsl_rng_uniform(rnd));
  }
  if(site[i] != site[h]){
    e2 = Exchange_Energy(site, neighbors, i, h);
    ie2 = Find_E_Index(state, e2);
    if(ie2 == -1){
      site[i] = -site[i];
      site[h] = -site[h];
      CurrentState = e2;
      CurrentIndex = StatesLength;
      if(CurrentState < MinState)
        MinState = CurrentState;	
      return(e2);
    }
  }

  if(site[i] == site[h]){
    state[CurrentIndex].logg += Logfnow;
		state[CurrentIndex].H += 1.0;	
  }
	else if(state[CurrentIndex].logg > state[ie2].logg || exp(state[CurrentIndex].logg - state[ie2].logg) > gsl_rng_uniform(rnd)){
		site[i] = -site[i];
    site[h] = -site[h];
		CurrentState = e2;
    CurrentIndex = ie2;
    if(CurrentState < MinState)
      MinState = CurrentState;	
		state[ie2].logg += Logfnow;
		state[ie2].H += 1.0;		
  }
  */
	return(-1);
}


/*	Check the flatness of the histogram	*/
int CheckFlat(struct StateInfo state[]){
	int i;
	double meanH, totalH=0.0, minH, maxH;
	minH = state[0].H;
	maxH = state[0].H;
 	totalH = maxH;
	for(i = 1; i<StatesLength; i++){
		totalH += state[i].H;
		if(state[i].H < minH)
      minH = state[i].H;
		else if(state[i].H > maxH) 
      maxH = state[i].H;
	}
	meanH = totalH / (float)StatesLength;
	if(minH > meanH*Flatness && meanH > maxH*Flatness){
    printf(" logf=%-.1e; Sweeps=%.1e; [minH,meanH,maxH]=[%-.1e, %-.1e, %-.1e]\n",  Logfnow/2.0, Total_Steps/(double)L, minH, meanH, maxH); 
		return(1);
	}
	if(Total_Steps > MaxSteps){
		Logfmin=1.0;
    printf("!!! MaxSteps %.2e reached! Logfmin is set to %.2f", MaxSteps, Logfmin);
	}
	return(0);
}


/* Reset the histogram to be 0*/
void Reset_H_0(struct StateInfo state[]){
    int i;
    for(i = 0; i < StatesLength; i++){
      state[i].H = 0.0;
    }
    return;
}

/* Swap two states: for the DIY quick sort*/
int Swap_States(struct StateInfo state[], int i, int j){
    struct StateInfo temp;
    temp.E = state[i].E;
    temp.logg = state[i].logg;
    temp.H = state[i].H;
    state[i].E = state[j].E;
    state[i].logg = state[j].logg;
    state[i].H = state[j].H;
    state[j].E = temp.E;
    state[j].logg = temp.logg;
    state[j].H = temp.H;
    return(0);
}
/* This function is just for quick sort*/
int Sort_Partition(struct StateInfo state[], int top, int bottom){
  int x = state[top].E;
  int i = top - 1;
  int j = bottom + 1;
  do{
      do{
          j--;
        }while(x > state[j].E );
      do{
          i++;
        }while(x < state[i].E);
      if(i < j) 
        Swap_States(state, i, j);
  }while(i < j);
  return(j);
}

/* Quick Sort the States by energy*/
void Quick_Sort_States(struct StateInfo state[], int top, int bottom){
    int middle;
    if(top < bottom){
        middle = Sort_Partition(state, top, bottom);
        Quick_Sort_States(state, top, middle);
        Quick_Sort_States(state, middle+1, bottom);
    }
    return;
}

/* Check the current energy state*/
int Check_CurrentEnergy(int sites[], int neighbors[], int neighloc[]){
  int i, j, neig_sum = 0, energy = 0;
  for(i = 0; i < (L+1); i++){
      neig_sum = 0;
      for(j = neighloc[i]; j < neighloc[i+1]; j++ ){
        neig_sum = neig_sum + sites[neighbors[j]];
      }
      energy = energy + sites[i] * neig_sum;
  }
  energy = energy / 2;
  return(energy);
}

/* Check the current mag state*/
int Check_CurrentMag(int sites[]){
  int i, mag = 0;
  for(i = 0; i < (L+1); i++){
    mag = mag + sites[i];
  }
  return(mag);
}

int Fit_Config(int site[], short configs[], short *config_num)
{
    int i=0, j=0, flag=0, diff_num=0;
    for(i = 0; i < *config_num; i++)
    {
        for (j = 0; j < (L+1); j++)
        {
            if(site[j] != configs[j+i*(L+1)])
            {
                  flag = 1;
                  break;
            }
        }
        if(flag==0)
          return 0;
        else
          diff_num++;
        flag=0;
    }
    if(*config_num==0 || diff_num == *config_num)
    {
        for(i=0; i< (L+1); i++)
        {
            configs[(*config_num)*(L+1) + i] = site[i];
        }
        (*config_num)++;
        printf("    config #%d for state=%d Added!\n", *config_num, CurrentState);
    }
    return 0;
}

int Check_Configs_Wrong(short configs[], short config_num)
{
    int ci=0, cj=0, i=0, flag=0;
    for(ci=0; ci<config_num-1; ci++)
    {
        for(cj=ci+1; cj<config_num; cj++)
        {
            for(i=0; i< (L+1); i++)
            {
                if(configs[ci*(L+1)+i] != configs[cj*(L+1)+i])
                {
                  flag = 1;
                  break;
                }
            }
            if(flag==1)
              flag=0;
            else
            {
              printf(" Row #%d == Row #%d (counting from 1)\n", ci+1, cj+1);
              return 1;
            }
        }
    }
    return 0;
  
}

/*	Function to read input parameters	*/
void Commandlineparse(int argc, char **argv, int *L, int *Flat,int *K, int *Seed, int *Name){
  	int i;
  	*Seed = getpid();
  	for (i = 1; i < argc; i++){  //Start at i = 1 to skip the command name.
    		if (argv[i][0] == '-'){
      			switch (argv[i][1]){
      				case 'L':       *L = atoi(argv[++i]);
        			break;
      				case 'f':       *Flat = atoi(argv[++i]);
        			break;
      				case 's':       *Seed = atoi(argv[++i]);
        			break;
      				case 'k':       *K = atoi(argv[++i]);
        			break;
      				case 'n':       *Name=atoi(argv[++i]);
        			break;
      				default:
        			fprintf(stderr,"\nError:  Incorrect option %s\n",argv[i]);
        			fprintf(stderr,"\nAvailable options: \n\
  				-L = L (HN3 length; it has to be 2^n, n=2,3,4,...)\n\
  				-k = highest level (HN3 length is 2^k))\n\
  				-f = flatness (Default = 80; or any integer in (0,99))\n\
  				-s = seed: random seed  (default=pid)\n\
    				\n");
        			exit(0);
      			}//switch
    		}//if
  	}//for
}



int main(int argc, char *argv[]){
  
  	int i, j, del1=0;	// General varialbe for loops
  	int new_state_energy;
  	char dos_file[24]="HPAFDOS_1D";	// Name of the output file of density of states
    char config_file[24]="HPAF_Config_1D";
  	char buf[8];
    struct StateInfo *state;
  	Commandlineparse(argc, argv, &L, &Flat, &K, &Seed, &Name);

  	rnd =  gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rnd, (unsigned long int)Seed);
    for(j = 0; j <100; j++) {
      		gsl_rng_uniform(rnd); // initial cycling on the random generator
  	}

  	printf("\n    ************ HN3 AF Ising Wang Landau Running ************ \n");
  	printf("\n      argc: %d\n      ",argc);
  	for(j = 0; j<argc; j++){
    		printf("%s ",argv[j]);
  	}
  
  	/* Setup the simulation according to the inputs*/
  	if(L==0){
    		L=Intpow(2,K);
  	}
  	else{
    		K=Intlog2(L);
  	}
  	Flatness=(float)Flat/100.0;
  	printf("\n\n L=%d; Flatness=%f;\n\n", L, Flatness);
  	MaxLength = (5*L); // this is just for HN3
  	int site[L+1];
    int neighbors[Intpow(2,K+2)+Intpow(2,K)], neighloc[L+2];
    short configs[32*(L+1)], config_num=0;
    for (i = 0; i < 32*(L+1); i++)
        configs[i] = -1;
    
    Initial_Steps= (double)((double)K*200*(double)L*(double)L);
    Check_Every_Steps=(double)((double)(K*K)*10.0*(double)L);
    state=(struct StateInfo *)malloc(MaxLength*sizeof(struct StateInfo));
    
    Initialize(state, site, neighbors, neighloc);
    printf(" Initialization Done \n");
  	Total_Steps=0.0;

  	sprintf(buf,"%d", K);
    strcat(dos_file, buf);
    strcat(config_file, buf);
    strcat(dos_file, "_r");
    strcat(config_file, "_r");
    sprintf(buf, "%d", Name);
    strcat(dos_file, buf);
    strcat(config_file, buf);
    
  	FILE *fp;

  	while(Total_Steps < Initial_Steps){
      Update(state, site, neighbors, neighloc);
      if(gsl_rng_uniform(rnd) < 0.01 && del1 < 3){
        if(Check_CurrentEnergy(site, neighbors, neighloc) == CurrentState)
          printf(" CurrentEnergy Correct = %d \n", CurrentState);
        else
          printf(" Error: CurrentEnergy is not right! \n CurrentState = %d,\
           actually %d \n", CurrentState, Check_CurrentEnergy(site, neighbors, neighloc));
        del1++;
      }
    }
    //printf(" Init Done. StatesLength = %d\n", StatesLength);
    //Print_States(state);
    //state = realloc(state, StatesLength*sizeof(struct StateInfo));
    if(state == NULL){
      fprintf(stderr, "\n Error: Memory reallocating of state(E,M,logg, H)!\n"); 
      return(2);
  	}
    printf(" Realloc Done\n");
    printf(" After Initial Steps %.2e, StatesLength = %d, MinEnergy = %d\n", Initial_Steps, StatesLength, MinState);
    Quick_Sort_States(state, 0, StatesLength-1);
    printf(" Quick_Sort Done \n");
    //Print_States(state);
    Inner_Loop_Steps=0.0;
    del1 = 0;
  	while(Logfnow > Logfmin){
      Inner_Loop_Steps += 1.0;
      new_state_energy = Update(state, site, neighbors, neighloc);
      if(gsl_rng_uniform(rnd) < 0.0001 && del1 < 3){//could delete
        if(Check_CurrentEnergy(site, neighbors, neighloc) == CurrentState)
          printf(" After sort, CurrentEnergy Correct = %d \n", CurrentState);
        else
          printf(" Error: CurrentEnergy is not right! \n CurrentState = %d,\
           actually %d \n", CurrentState, Check_CurrentEnergy(site, neighbors, neighloc));
        del1++;
      }
      if(new_state_energy != -1){	// if new state found, we add another states
        StatesLength++;
        if(StatesLength > MaxLength) {
          Print_States(state);
          fprintf(stderr, "\n Error: the states length %d is too long;\n Check the States above!!!\n", StatesLength);
          exit(0);
        }
        printf(" New State: %d: coming; StatesLength %d\n", new_state_energy, StatesLength);
        /* These are to resize the matrix size; which is not necessary. just waste at most several Mb memory.
        state = (struct StateInfo *)realloc(state, StatesLength*sizeof(struct StateInfo));
        if(state == NULL){
            fprintf(stderr, "\n Error: Memory reallocating of state(E, logg, H)!\n"); 
            return(2);
        } */
        for(i = StatesLength-2; i > 0; i--){
            if(state[i].E > new_state_energy){
                state[i+1].E = new_state_energy;
                state[i+1].logg = state[i].logg * 0.8;
                state[i+1].H = 0.0;
                Logfnow = Logfnow * 8;
                if(Logfnow > 1) Logfnow = 1;
                Reset_H_0(state);
                break;
            }
            else{
                state[i+1].E = state[i].E;
                state[i+1].logg = state[i].logg;
                state[i+1].H = state[i].H;
            }
        }
        printf(" New State: %d added to list; Min State: %d; CurrentState = %d\n", new_state_energy, MinState, CurrentState);
        //Print_States(state);
        printf(" \n\n");
        Inner_Loop_Steps=0.0;
        
      }
      if(Inner_Loop_Steps > Check_Every_Steps){
        Inner_Loop_Steps = 0.0;
        if(CheckFlat(state)){
          Reset_H_0(state);
          Logfnow=Logfnow/2.0;
        }
      }
      
      
      // check configurations if the lowest states is likely found
      if(new_state_energy != -1 && Logfnow < 0.125)
      {
          for(i=0; i < config_num*(L+1); i++)
            configs[i] = -1;
          config_num = 0;
      }
      else if(Logfnow < 0.125 && CurrentState == MinState && config_num<32)
          Fit_Config(site, configs, &config_num);
  	}
 
    printf("\n MinState: %d; \n MaxState: %d; \n States Length: %d; \n Total MC steps: %-.4e;\
     \n Total random steps: %-.4e; \n ", 
     MinState, MaxState, StatesLength, Total_Steps/(double)(L), Total_Steps);
  	fp=fopen(dos_file,"w");
  	MinStateLogg=state[0].logg;
  	for(i = 0; i < StatesLength; i++){
      state[i].logg = state[i].logg-MinStateLogg + log(2.0);
      fprintf(fp, "%-8d %-18.8f\n", state[i].E, state[i].logg);
  	}
  	fclose(fp);
    
    if(Check_Configs_Wrong(configs, config_num))
    {
        printf("Error: Repeated Config Found!\n");
    }
    strcat(config_file, "_c");
    sprintf(buf, "%d", config_num);
    strcat(config_file, buf);
    fp=fopen(config_file, "w");
    fprintf(fp, " HNNP\n L : %d;\n MinState: %d; \n MinStates Found: %d; \n MaxState: %d; \n States Length: %d; \n Total MC steps: %-.4e;\
     \n Total random steps: %-.4e; \n\n", 
     L, MinState, config_num, MaxState, StatesLength, Total_Steps/(double)(L), Total_Steps);
    fprintf(fp, "*********** HNNP, L=%d, Lowest Energy=%d, All configurations *************\n", L, MinState);
    for(j=0; j< config_num; j++)
    {
        if(j==0)
           fprintf(fp, "%-10s", "site #");
        fprintf(fp, "C#%-6d  ", j+1);
    }
    fprintf(fp, "\n");
    for(i=0; i< (L+1); i++)
    {
      for (j = 0; j < config_num; j++)
      {
          if(j==0)
            fprintf(fp, "%-9d ", i);
          fprintf(fp, "%-+9d ", configs[j*L + i]);
      }
      fprintf(fp, "\n");
    }
      
    fclose(fp);
    
  /* code below are Calculations: 
   * we skip this part and do this using Python 
   * so that we can visuluize everything 
   */
	/*
	for(i=0; i<INV_MU_L; i++){
		z[i]=0.0;
		rho[i]=0.0;
		s[i]=0.0;
		shannon[i]=0.0;
	}
	*/
	/*	Calculating Z, packing fraction, entropy and so on */
	/*
	for(i=0; i<INV_MU_L; i++){
		for(j=0; j<StatesLength; j++){
			z_temp = expl((long double)state[j].logg)*expl((long double)(1/InvMu[i]*(long double)j));
			z[i] += z_temp;
			rho[i] += (long double)(j*z_temp);
		}
		for(j=0; j<StatesLength; j++){
			z_temp= expl((long double)state[j].logg)*expl((long double)(1/InvMu[i]*(double)j));
			pj=(long double)(z_temp/z[i]);
			if(pj>1e-12){
				shannon[i] -= pj*log2l(pj);
			}
		}
		rho[i]=(long double)(rho[i]/z[i]/(double)L);
		s[i]=1.0/(long double)L*logl(z[i])-1/InvMu[i]*rho[i];
	}

  	fp=fopen(pes_file,"w");
	for(i=0; i<INV_MU_L; i++){
		fprintf(fp, "%-8.4f  %-18.9LG %-18.9LG  %-18.9LG \n", InvMu[i], rho[i], s[i], shannon[i]);
	}
	fclose(fp);
	*/
 	return(0);
}
