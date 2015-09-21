#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define PI 3.14159265
#define TWOPI 6.283185307
int main(int argc, char **argv)
{
int a;
double **histogram;
double radial;
double height;
double phi_slice;
double phitolerance;
double phitolOVERtwo;
double halfheight;
long i;
long j;
long r_num;
long r_bin;
long z_num;
long z_bin;
double r_cyn;
double z;
double dr;
double dz;
double r;
double th;
double ph;
FILE *file;
FILE *data;
FILE *output;
char filename[18];
i = 0;
while(i < argc){
    if ( strcmp(argv[i],
         radial = atof(argv[++i]);
    else if ( strcmp(argv[i],       "--h") == 0 )
         height = atof(argv[++i]);
    else if ( strcmp(argv[i],    "--ph") == 0 )
         phi_slice = atof(argv[++i]);
    else if ( strcmp(argv[i],       "--dr") == 0 )
         dr = atof(argv[++i]);
    else if ( strcmp(argv[i],       "--dz") == 0 )
         dz = atof(argv[++i]);
         else if ( strcmp(argv[i],       "--dph") == 0 )
        phitolerance = atof(argv[++i]);
   else if ( strcmp(argv[i], "--help") == 0 ){
        fprintf(stderr, "\nStrong Field Cross Section Plotter\n"
                        "V1.0  Author: Carl Rodriguez\n\n"
                        "Calculates an azimuthal subsection of an accretion disk\n"
"Designed for the near field regime, close to the central body\n\n"
        exit(0);}
    i++;
}
"--r Radius of interest\n"
"--h Height of interest\n"
"--ph Location of phi slice (0 to pi)\n"
"--dr Radial bin size\n"
"--dz Height bin size\n"
"--dph Thickness of phi slice\n");
halfheight = height/2.;
phitolOVERtwo = phitolerance/2.;
r_num = (long)(ceil((2*radial)/dr));
z_num = (long)(ceil(height/dz));
histogram = (double **)malloc(z_num * sizeof(double *));
for( i = 0 ; i < z_num ; i++)
    histogram[i] = (double *)malloc(r_num * sizeof(double));
for( i = 0 ; i < z_num ; i++ )
for( j = 0 ; j < r_num ; j++ )
   histogram[i][j] = 0;
file = fopen("interesting.full.dat","r");
if (file==0) {
    fprintf(stderr, "Cannot open interesting file list, aborting\n");
return 0; }
while (!feof(file)) {
    fscanf(file, "%d\n", &a);
    sprintf(filename, "Kerrpolar_%d.dat", a);
    data = fopen(filename, "r");
    if (data==0) {
        fprintf(stderr, "Cannot open %s, skipping\n", filename);
    }
    else{
        while (!feof(data)) {
            fscanf(data, "%lf       %lf     %lf\n", &r, &th, &ph);
            if (fmod(fabs(ph-phitolOVERtwo+phi_slice),TWOPI) < phitolerance){
                r_cyn = radial + (r * sin(th));
                z = r * cos(th) + halfheight;
                r_bin = (long)floor(r_cyn / dr);
                z_bin = (long)floor(z / dz);
                if((fabs(r_bin) < r_num) && (z_bin < z_num) && (r_bin >= 0 )
&& (z_bin >= 0))
} }
    }
fclose(data);
}
fclose(file);
        histogram[z_bin][r_bin] += 1;
}
else if(fmod(fabs(ph-phitolOVERtwo+phi_slice) + PI,TWOPI)
  < phitolerance){
    r_cyn = radial - (r * sin(th));
    z = r * cos(th) + halfheight;
    r_bin = (long)floor(r_cyn / dr);
    z_bin = (long)floor(z / dz);
    if((fabs(r_bin) < r_num) && (z_bin < z_num) && (r_bin >= 0 ) && (z_bin >= 0))
    histogram[z_bin][r_bin] += 1;
    output = fopen("StrongHistogram.dat", "w");
    for(i = 0 ; i < z_num ; i++){
        for(j = 0 ; j < r_num ; j++)
            fprintf(output, "%lf\t", histogram[i][j]);
            fprintf(output, "\n");
}
fclose(output);
return 0;
}