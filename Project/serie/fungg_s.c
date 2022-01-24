/*
   CA - OpenMP
   fungg_s.c
   Routines used in gengroups_s.c program

   TO BE COMPLETED
***************************************************************/


#include <math.h>
#include <float.h>
#include "definegg.h"		// definition of constants


/* 1 - Function to calculate the genetic distance; Euclidean distance between two elements.
       Input:   two elements of NFEAT characteristics (by reference)
       Output:  distance (double)
***************************************************************************************************/

double geneticdistance (float *elem1, float *elem2)
{
   double cum = 0;
   for(int i =0; i<NFEAT; i++){

      double ans = elem1[i] - elem2[i];
      cum = cum + pow(ans, 2);
   }
   return sqrt(cum);

}



/* 2 - Function to calculate the closest group (closest centroid) for each element.
   Input:   nelems   number of elements, int
            elems    matrix, with the information of the elements, of size MAXELE x NFEAT, by reference
            cent    matrix, with the centroids, of size NGROUPS x NFEAT, by reference
   Output:  grind   vector of size MAXELE, by reference, closest group for each element
***************************************************************************************************/

void closestgroup (int nelems, float elems[][NFEAT], float cent[][NFEAT], int *grind)
{
   int g;
   double tdis, mind;
   for(int i = 0; i<nelems;i++){ 
      mind = DBL_MAX;
      for(int j=0; j<NGROUPS;j++){
         tdis = geneticdistance(elems[i],cent[j]);
         if(tdis < mind){
            mind = tdis;
            g = j;
         }
      }
      grind[i] = g;
   }
}




/* 3 - Function to calculate the compactness of each group (average distance between all the elements in the group) 
   Input:  elems     elements (matrix of size MAXELE x NFEAT, by reference)
           iingrs   indices of the elements in each group (vector of size NGROUPS with information for each group)
   Output: compact  compactness of each group (vector of size NGROUPS, by reference) 
***************************************************************************************************/

void groupcompactness (float elems[][NFEAT], struct ginfo *iingrs, float *compact)
{
   int o, a;
   double cum, cont;

   for(int i =0; i < NGROUPS; i++){
      int telem = iingrs[i].size;
      if(telem < 2){
         compact[i] = 0;
      }else{
         cum = 0.0;
         cont = 0.0;
         for(int j = 0; j < telem; j++){
            a = iingrs[i].members[j];
            for( int k = j+1; k<telem; k++){
               o = iingrs[i].members[k];
               cum = cum + geneticdistance(elems[a], elems[o]);
               cont = cont + 1.0;
            }
         }
         compact[i] = (float) (cum/cont);
      }
   }
}




/* 4 - Function to analyse diseases 
   Input:  iingrs   indices of the elements in each group (matrix of size NGROUPS x MAXELE, by reference)
           dise     information about the diseases (NGROUPS x TDISEASE)
   Output: disepro  analysis of the diseases: maximum, minimum of the medians and groups
***************************************************************************************************/

void diseases (struct ginfo *iingrs, float dise[][TDISEASE], struct analysis *disepro)
{
   int a, j, k, gmax, telem, gmin;
   float bmedian, lmedian, tmedian;
   for(int i = 0; i < TDISEASE; i++){
      lmedian = FLT_MAX;
      bmedian = FLT_MIN;
      for(j = 0; j < NGROUPS; j++){
         telem = iingrs[j].size;
         float pro[telem];
         for( k = 0; k<telem; k++){ 
            a = iingrs[j].members[k];
            pro[k] = dise[a][i];
         }
         quicksort(&pro, 0, telem-1);   

         int size = sizeof pro/sizeof pro[0];
         if(size>1){ 
         tmedian = pro[(telem/2)];
         }else if(size ==1){
            tmedian = pro[0];
         }
         if(tmedian > 0 && tmedian < lmedian){
            lmedian = tmedian;
            gmin = j;
         } else if(tmedian > bmedian){
               bmedian = tmedian;
               gmax = j;
         }
      }
      tmedian = 0.0;
      disepro[i].gmax = gmax;
      disepro[i].mmin = lmedian;
      disepro[i].gmin = gmin;
      disepro[i].mmax = bmedian;
   }
}


// Swapping algorithm
static inline
void swap(float *a, float *b)
{
    float temp = *a;
    *a = *b;
    *b = temp;
}


// Partitioning algorithm
static
int partition(float *L, int left, int right)
{
    int pivot = left;
      float p_val = L[pivot];

    while (left < right)
    {
        while (L[left] <= p_val)
            left++;
        while (L[right] > p_val)
            right--;
        if (left < right)
            swap(&L[left], &L[right]);
    }
    swap(&L[pivot], &L[right]);
    return right;
}

// Quicksort recursion

void quicksort(float *L, int start, int end)
{
    if (start >= end)
        return;
    //dump_list("PRE-PARTITION", L, start, end);
    int splitPoint = partition(L, start, end);
    //dump_list("POST-PARTITION", L, start, end);
    //printf("Split point: %d\n", splitPoint);
    quicksort(L, start, splitPoint - 1);
    quicksort(L, splitPoint + 1, end);
}