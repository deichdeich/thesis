#include "RKheader.h"

/*Converts from the relativistic constants of motion (energy, spin, and carter’s constant)
 to the starting coordinates for the simulation*/

int convert
(
  double energy,
  double jz,
  double carter,
  double M,
  double a,
  double *f
){
 double test;
  f[1] = -energy - (2*M*f[2]*(a*(a*energy + jz) + energy*pow(f[2],2)))/((pow(a,2)*pow(cos(f[4]),2) + pow(f[2],2))*(pow(a,2) - 2*M*f[2] + pow(f[2],2)));
  
  f[5] = sqrt(2*a*energy*jz + carter - pow(jz,2)*pow((1/sin(f[4])),2) + pow(a,2)*(pow(cos(f[4]),2) - pow(energy,2)*pow(sin(f[4]),2)))/ pow(pow(f[2],2) + pow(a,2)*pow(cos(f[4]),2),2);
  
  test = 2*a*energy*jz + carter - pow(jz,2)*pow((1/sin(f[4])),2) + pow(a,2)*(pow(cos(f[4]),2) - pow(energy,2)*pow(sin(f[4]),2));
  
  if(test <= 0.)f[5] = 0.;
    /*shody piece of chicanery; it handles floating point errors that somtimes
      conspire to give an imaginary theta_dot for planar orbits*/
  f[7] = (pow(a,2)*jz*pow((1/tan(f[4])),2) + f[2]*(-2*a*energy*M + jz*pow((1/sin(f[4])),2)*(-2*M + f[2])))/((pow(a,2)*pow(cos(f[4]),2) + pow(f[2],2))*(pow(a,2) - 2*M*f[2] + pow(f[2],2)));
  
  return 0;
}
/*The first order vector form of the Kerr Lagrangian (divided by tdot to parameterize in
coordinate time)*/
int func
(
  double t,
  double tdot,
  double r,
  double rdot,
  double th,
  double thdot,
  double ph,
  double phdot,
  double *g,
  double M,
  double a,
  double b,
  double disk_bound
){
  g[0] = 1;
  g[1] = (M*rdot*(2.*(pow(a,2.) + pow(r,2.))*tdot*(-pow(r,2.) + pow(a,2.)*pow(
         cos(th),2.)) + a*phdot*(6.*pow(r,4.) - 2.*pow(a,4.)*pow(cos(th),2.) +
         pow(a,2.)*pow(r,2.)*(3. + cos(2.*th)))*pow(sin(th),2.)) - 2.*pow(a,2.)
         *M*r*(pow(a,2.) - 2.*M*r + pow(r,2.))*thdot*(-tdot + a*phdot*pow(sin(
         th),2.))*sin(2.*th))/((pow(a,2.) - 2.*M*r + pow(r,2.))*tdot*pow(pow(r,
         2.) + pow(a,2.)*pow(cos(th),2.),2.));
  g[2] = rdot/tdot;
  g[3] = -((pow(a,2.) - 2.*M*r + pow(r,2.))*((-2.*r*pow(rdot,2.))/(pow(a,2.) -
￼2.*M*r + pow(r,2.)) - 2.*r*pow(thdot,2.) + (2.*M*pow(tdot,2.)*(pow(r,2.) - pow(a,2.)*pow(cos(th),2.))) / pow(pow(r,2.) + pow(a,2.)*pow(cos(th),
     2.),2.) + (2.*(M - r)*pow(rdot,2.)*(pow(r,2.) + pow(a,2.)*pow(cos(th),2.
     ))) / pow(pow(a,2.) - 2.*M*r + pow(r,2.),2.) - (8.*a*M*phdot*pow(r,2.)*
     tdot*pow(sin(th),2.)) / pow(pow(r,2.) + pow(a,2.)*pow(cos(th),2.),2.) +
     (4.*a*M*phdot*tdot*pow(sin(th),2.))/(pow(r,2.) + pow(a,2.)*pow(cos(th),2.
     )) - (pow(phdot,2.)*pow(sin(th),2.)*(4.*pow(r,3.) + pow(a,2.)*r*(3. +
     cos(2.*th)) + 2.*pow(a,2.)*M*pow(sin(th),2.)))/(pow(r,2.) + pow(a,2.)*
     pow(cos(th),2.)) + (2.*pow(phdot,2.)*r*pow(sin(th),2.)*(pow(pow(a,2.) +
     pow(r,2.),2.) - pow(a,2.)*(pow(a,2.) - 2.*M*r + pow(r,2.))*pow(sin(th),2.
     ))) / pow(pow(r,2.) + pow(a,2.)*pow(cos(th),2.),2.) -(2*rdot*(-2*r*rdot +
     pow(a,2)*thdot*sin(2*th)))/(pow(a,2) - 2*M*r + pow(r,2)))/(2.*tdot*(pow(
     r,2.)+ pow(a,2.)*pow(cos(th),2.))) - (1./(1.+exp(1.*(r - disk_bound))))
     *b*rdot/r;
g[4] = thdot/tdot;
g[5] = -((2.*r*rdot*thdot + ((-(pow(phdot,2.)*pow(pow(a,2.) - 2.*M*r + pow(r,
       2.),2.)) + pow(a,2.)*(pow(rdot,2.) - (pow(a,2.) - 2.*M*r + pow(r,2.))*
       pow(thdot,2.)))*sin(2.*th)) / (2.*(pow(a,2.) - 2.*M*r + pow(r,2.))) -
       (M*r*pow(phdot*(pow(a,2.) + pow(r,2.)) - a*tdot,2.)*sin(2.*th))/pow(pow
       (r,2.) + pow(a,2.)*pow(cos(th),2.),2.))/(tdot*(pow(r,2.) + pow(a,2.)*
       pow(cos(th),2.)))) - (1./(1.+exp(1.*(r - disk_bound))))*b*thdot/r;
g[6] = phdot/tdot;
g[7] = -(csc(th)*(-4.*a*M*tdot*(8.*r*(pow(a,2.) - 2.*M*r + pow(r,2.))*thdot*
       cos(th) + 4.*rdot*(-pow(r,2.) + pow(a,2.)*pow(cos(th),2.))*sin(th)) +
       phdot*(16.*(pow(a,2) - 2.*M*r + pow(r,2.))*thdot*cos(th)*(pow(pow(r,2.)
       + pow(a,2.)*pow(cos(th),2.),2.) + 2.*pow(a,2.)*M*r*pow(sin(th),2.)) +
       4.*rdot*sin(th)*(-8.*M*pow(r,4.) + 4.*pow(r,5.) + 8.*pow(a,2.)*pow(r,3.)
       *pow(cos(th),2.) + 4.*pow(a,4.)*r*pow(cos(th),4.) - 2.*pow(a,2.)*M*pow
       (r,2.)*(3. + cos(2.*th)) + pow(a,4.)*M*pow(sin(2.*th),2.))))) / (8.*(
       pow(a,2.) - 2.*M*r + pow(r,2.))*tdot*pow(pow(r,2.) + pow(a,2.)*pow(
       cos(th),2.),2.)) - (1./(1.+exp(1.*(r - disk_bound))))*b*phdot/r;
return 0; }
double kerr_error
(
  double *f_temp,
  double *error,
  double M,
  double a
){
  double proper_error;
  double proper_v_error;
  proper_error = (-1 + (2*M*f_temp[2])/(pow(f_temp[2],2) + pow(a,2)*pow(cos(
         f_temp[4]),2)))*pow(error[0],2) + ((pow(f_temp[2],2) + pow(a,2)*pow(cos(
         f_temp[4]),2))/(pow(a,2) - 2*M*f_temp[2] + pow(f_temp[2],2)))*pow(error[2],
         2) + (pow(f_temp[2],2) + pow(a,2)*pow(cos(f_temp[4]),2))*pow(error[4], 2) +
         ((pow(sin(f_temp[4]),2)*(pow(pow(a,2) + pow(f_temp[2],2),2) - pow(a,2)*(
         pow(a,2) - 2*M*f_temp[2] + pow(f_temp[2],2))*pow(sin(f_temp[4]),2)))/(pow(
         f_temp[2],2) + pow(a,2)*pow(cos(f_temp[4]),2)))*pow(error[6], 2) + ((-2*a*M
         *f_temp[2]*pow(sin(f_temp[4]),2))/(pow(f_temp[2],2) +  pow(a,2)*pow(cos(
         f_temp[4]),2)))*error[0]*error[6];
          proper_v_error = (-1 + (2*M*f_temp[2])/(pow(f_temp[2],2) + pow(a,2)*pow(cos(
         f_temp[4]),2)))*pow(error[1],2) + ((pow(f_temp[2],2) + pow(a,2)*pow(cos(
         f_temp[4]),2))/(pow(a,2) - 2*M*f_temp[2] + pow(f_temp[2],2)))*pow(error[3],
         2) + (pow(f_temp[2],2) + pow(a,2)*pow(cos(f_temp[4]),2))*pow(error[5], 2) +
         ((pow(sin(f_temp[4]),2)*(pow(pow(a,2) + pow(f_temp[2],2),2) - pow(a,2)*(
         pow(a,2) - 2*M*f_temp[2] + pow(f_temp[2],2))*pow(sin(f_temp[4]),2)))/(pow(
         f_temp[2],2) + pow(a,2)*pow(cos(f_temp[4]),2)))*pow(error[7], 2) + ((-2*a*M
         *f_temp[2]*pow(sin(f_temp[4]),2))/(pow(f_temp[2],2) + pow(a,2)*pow(cos(
         f_temp[4]),2)))*error[1]*error[7];
  return sqrt(pow(proper_error, 2.) + pow(proper_v_error, 2.));
}


int main(int argc, char **argv)
{
  double M;
  double a;
  double b;
  double disk_bound;
  double energy;
  double energy_temp;
  double jz;
  double phdot_central;
  double carter;
  double carter_min;
  double *f;
  double *g;
  double *f_temp;
  double *error;
  double ds;
  double ds_old;
  double t;
  double dt;
  double input_error;
  double proper_error;
  double proper_v_error;
  double total_error;
  double event_horizon;
  double horizon_tolerance;
  double t_max;
  double k1[8], k2[8], k3[8], k4[8], k5[8], k6[8];
  int i, j;
  long NN;
  FILE *filefull, *filepolar;
  char Kerrfull[20], Kerrpolar[20];
  /*Sets up the GSL random number generator*/
  const gsl_rng_type * T;
  gsl_rng * r;
  T = gsl_rng_default;
  gsl_rng_env_setup();
  r = gsl_rng_alloc (T);
/*handles the in/out vectors as pointers for purposes of function passing*/
g = (double *)malloc(8.*sizeof(double));
f = (double *)malloc(8.*sizeof(double));
f_temp = (double *)malloc(8.*sizeof(double));
error = (double *)malloc(8.*sizeof(double));
for( i = 0 ; i < 8 ; i++){
  g[i] = 0.;
  f[i] = 0.;
  f_temp[i] = 0.;
  error[i] = 0.;
  k1[i] = 0.;
  k2[i] = 0.;
  k3[i] = 0.;
  k4[i] = 0.;
  k5[i] = 0.;
  k6[i] = 0.;
}
/*defauly error and event horizon tolerance*/
input_error = 0.00001;
horizon_tolerance = 0.01;
dt = 1.;
t = 0.;
i = 0;
j = 0;
f[4] = PIover2;
NN = 0;
b = 0.;
/*The input/output for the command line*/
while(i < argc){
  if ( strcmp(argv[i],  "--M") == 0 )
    M = atof(argv[++i]);
  else if ( strcmp(argv[i],  "--a") == 0 )
    a = atof(argv[++i]);
  else if ( strcmp(argv[i],  "--b") == 0 )
    b = atof(argv[++i]);
  else if ( strcmp(argv[i],  "--db") == 0 )
    disk_bound = atof(argv[++i]);
  else if ( strcmp(argv[i],  "--r") == 0 )
    f[2] = atof(argv[++i]);
  else if ( strcmp(argv[i],  "--e") == 0 )
    energy = atof(argv[++i]);
  else if ( strcmp(argv[i],  "--jz") == 0 )
    jz = atof(argv[++i]);
  else if ( strcmp(argv[i],  "--k") == 0 )
    carter = atof(argv[++i]);
  else if ( strcmp(argv[i],  "--th") == 0 )
    f[4] = atof(argv[++i]);
  else if ( strcmp(argv[i],  "--ph") == 0 )
    f[6] = atof(argv[++i]);
  else if ( strcmp(argv[i],  "--ie") == 0 )
    input_error = atof(argv[++i]);
  else if ( strcmp(argv[i],  "--eh") == 0 )
   horizon_tolerance = atof(argv[++i]);
  else if ( strcmp(argv[i],  "--tm") == 0 )
    t_max = atof(argv[++i]);
  else if ( strcmp(argv[i],  "--dt") == 0 )
    dt = atof(argv[++i]);
  else if ( strcmp(argv[i],  "--nn") == 0 )
    NN = atoi(argv[++i]);
  else if ( strcmp(argv[i],  "--help") == 0 ){
    fprintf(stderr, "\nKerr Geodesic Simulator V2.0\n"
      "Author: Carl Rodriguez\n\n"
      "Uses the relativistic lagrangian and the full kerr metric \n"
      "to compute the trajectories of a test body with an RK54 ODE solver\n\n"
      "Usage: --M  mass of central body\n"
      "   --a   spin of central body\n"
      "   --b   damping coefficient (default is 0)\n"
      "   --db  disk boundary\n"
      "   --r   initial radial position\n"
      "   --e   orbital energy\n"
      "   --jz  orbital spin\n"
      "   --k   Carter’s constant\n"
      "   --th  inital theta position (default is pi/2)\n"
      "   --ie  Error epsilon (default is 1e-5)\n"
      "   --ph  initial phi position (default is 0)\n"
      "   --eh  Event Horizon epsilon (default is 0.01)\n"
      "   --dt  Time Intervals (default is 1m)\n"
      "   --tm  Max run time (in coordinate time)\n"
      "   --nn  Run/RNG seed number (for cluster)\n");
      exit(0);}
i++; }
sprintf(Kerrfull, "Kerrfull_%ld.dat", NN);
sprintf(Kerrpolar, "Kerrpolar_%ld.dat", NN);
gsl_rng_set (r, NN);
/*Adds gaussian variation to the input parameters, then converts them to
 spherical initial conditions. Also checks to make sure Carter’s Constant
 is physical*/
if(NN != 0){
  phdot_central = (4*(pow(f[2],2) + pow(a,2)*pow(cos(f[4]),2))*pow(csc(f[4]),2)*
    (jz*f[2]*(-2*M + f[2]) + a*(a*jz*pow(cos(f[4]),2) + 2*energy*M*f[2]*pow(sin(
    f[4]),2)))) / ((pow(a,2) + f[2]*(-2*M + f[2]))*pow(pow(a,2) + 2*pow(f[2],2)
    + pow(a,2)*cos(2*f[4]),2));
  do energy_temp = gsl_ran_gaussian(r, ENERGY_SIGMA);
    while(energy_temp + energy <= 0.);
  energy += energy_temp;
  f[4] += gsl_ran_gaussian(r, THETA_SIGMA);
  f[2] += gsl_ran_gaussian(r, R_SIGMA);
  jz = ((-4*a*f[2]*energy*M + (pow(a,2) + 2*pow(f[2],2))*(pow(a,2) + f[2]*(f[2] -
       2*M))*phdot_central + pow(a,2)*(pow(a,2) + f[2]*(f[2] -2*M))*phdot_central
       *cos(2*f[4]))*pow(sin(f[4]),2))/(pow(a,2) + 2*f[2]*(f[2] - 2*M) + pow(a,2)
       *cos(2*f[4]));
       
  carter += gsl_ran_gaussian(r, CARTER_SIGMA);
}
fprintf(stderr, "RNG Seed:%ld R:%lf Jz:%lf e:%lf theta:%lf K:%lf\n",NN,
f[2],jz,energy,f[4],carter);
carter_min = -2.*a*energy*jz - pow(a,2.)*pow(cos(f[4]),2.) + pow(jz,2.)*pow(csc(
       f[4]),2.) + pow(a,2.)*pow(energy,2.)*pow(sin(f[4]),2.);
if(carter < carter_min){
    fprintf(stderr,"WARNGING: Input Carter’s Constant is below possible minimum.\n"
                  "Setting constant to %0.15lf\n", carter_min);
    carter = carter_min;
}
convert(energy, jz, carter, M, a, f);
filefull = fopen(Kerrfull, "w");
filepolar = fopen(Kerrpolar, "w");
ds = input_error;
event_horizon = M + sqrt(pow(M,2.) - pow(a,2.));
i = 0;
/*RK54 ODE Solver*/
while(f[0] <= t_max){
  while(f[0] >= t + dt ){
    t += dt;
    energy = -(pow(a,2)*pow(cos(f[4]),2)*f[1] + f[2]*(f[1]*(-2*M + f[2]) + 2*a*M
             *f[7]*pow(sin(f[4]),2)))/(pow(a,2)*pow(cos(f[4]),2) + pow(f[2],2));
    jz = (pow(sin(f[4]),2)*(-2*a*M*f[1]*f[2] + pow(pow(a,2) + pow(f[2],2),2)*f[7]
         - pow(a,2)*(pow(a,2) + f[2]*(-2*M + f[2]))*f[7]*pow(sin(f[4]),2)))/(pow(
         a,2)*pow(cos(f[4]),2) + pow(f[2],2));
    fprintf(filefull, "%lf     %lf     %lf     %lf %0.12lf %lf\n",f[0], f[2],
         f[4], f[6], energy, jz);
}
  fprintf(filepolar, "%lf %lf %lf\n", f[2], f[4], f[6]);
  if( i == 0 && f[2] < disk_bound){
    fprintf(stderr, "Disk Geodesic!\n");
    fprintf(stdout, "%ld\n", NN);
    ++i;
}
  do{
    func(f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],g,M,a,b,disk_bound);
    for(j = 0 ; j < 8 ; j++) k1[j] = ds * g[j];
    func(f[0]+B_21*k1[0],f[1]+B_21*k1[1],f[2]+B_21*k1[2],f[3]+B_21*k1[3],f[4]+B_21
        *k1[4],f[5]+B_21*k1[5],f[6]+B_21*k1[6],f[7]+B_21*k1[7],g,M,a,b,disk_bound);
    for(j = 0 ; j < 8 ; j++) k2[j] = ds * g[j];
    func(f[0]+B_31*k1[0]+B_32*k2[0],f[1]+B_31*k1[1]+B_32*k2[1],f[2]+B_31*k1[2]+
        B_32*k2[2],f[3]+B_31*k1[3]+B_32*k2[3],f[4]+B_31*k1[4]+B_32*k2[4],f[5]+
        B_31*k1[5]+B_32*k2[5],f[6]+B_31*k1[6]+B_32*k2[6],f[7]+B_31*k1[7]+B_32*k2[7]
        ,g,M,a,b,disk_bound);
    for(j = 0 ; j < 8 ; j++) k3[j] = ds * g[j];
    func(f[0]+B_41*k1[0]+B_42*k2[0]+B_43*k3[0],f[1]+B_41*k1[1]+B_42*k2[1]+B_43*k3[1],f[2]+B_41*k1[2]+B_42*k2[2]+B_43*k3[2],f[3]+B_41*k1[3]+B_42*k2[3]+B_43*k3[3],
          f[4]+B_41*k1[4]+B_42*k2[4]+B_43*k3[4],f[5]+B_41*k1[5]+B_42*k2[5]+B_43*k3[5],
          f[6]+B_41*k1[6]+B_42*k2[6]+B_43*k3[6],f[7]+B_41*k1[7]+B_42*k2[7]+B_43*k3[7],
          g,M,a,b,disk_bound);
      for(j = 0 ; j < 8 ; j++) k4[j] = ds * g[j];
      func(f[0]+B_51*k1[0]+B_52*k2[0]+B_53*k3[0]+B_54*k4[0],f[1]+B_51*k1[1]+B_52*k2[1]
          +B_53*k3[1]+B_54*k4[1],f[2]+B_51*k1[2]+B_52*k2[2]+B_53*k3[2]+B_54*k4[2],f[3]
          +B_51*k1[3]+B_52*k2[3]+B_53*k3[3]+B_54*k4[3],f[4]+B_51*k1[4]+B_52*k2[4]+B_53
          *k3[4]+B_54*k4[4],f[5]+B_51*k1[5]+B_52*k2[5]+B_53*k3[5]+B_54*k4[5],f[6]+B_51*
          k1[6]+B_52*k2[6]+B_53*k3[6]+B_54*k4[6],f[7]+B_51*k1[7]+B_52*k2[7]+B_53*k3[7]
          +B_54*k4[7],g,M,a,b,disk_bound);
      for(j = 0 ; j < 8 ; j++) k5[j] = ds * g[j];
      func(f[0]+B_61*k1[0]+B_62*k2[0]+B_63*k3[0]+B_64*k4[0]+B_65*k5[0],f[1]+B_61*k1[1]+
          B_62*k2[1] + B_63*k3[1]+B_64*k4[1]+B_65*k5[1],f[2]+B_61*k1[2]+B_62*k2[2]+B_63*
          k3[2]+B_64*k4[2]+B_65*k5[2],f[3]+B_61*k1[3]+B_62*k2[3]+B_63*k3[3]+B_64*k4[3]+
          B_65*k5[3],f[4]+B_61*k1[4]+B_62*k2[4]+ B_63*k3[4]+B_64*k4[4]+B_65*k5[4],f[5]+
          B_61*k1[5]+B_62*k2[5]+B_63*k3[5]+B_64*k4[5]+B_65*k5[5],f[6]+B_61*k1[6]+B_62*
          k2[6]+B_63*k3[6]+B_64*k4[6]+B_65*k5[6],f[7]+B_61*k1[7]+B_62*k2[7]+B_63*k3[7]+
          B_64*k4[7]+B_65*k5[7],g,M,a,b,disk_bound);
      for(j = 0 ; j < 8 ; j++){
        k6[j] = ds * g[j];
        error[j] = CDIFF_1*k1[j] + CDIFF_3*k3[j] + CDIFF_4*k4[j] +
        CDIFF_5*k5[j] + CDIFF_6*k6[j];
        f_temp[j] = f[j] + k1[j]*C_1 + k3[j]*C_3 + k4[j]*C_4 + k6[j]*C_6;
}
      total_error = kerr_error(f_temp, error, M, a);
      ds_old = ds;
      ds *= pow(input_error / total_error, 0.2);
    }while(ds_old > ds);
    for(j = 0 ; j < 8 ; j++) f[j] = f_temp[j];
    if(f[2] - event_horizon < horizon_tolerance ){
      fprintf(stderr, "INSPIRAL!\n");
      break;
} }
  fclose(filefull);
  fclose(filepolar);
  if(i==0){
    remove(Kerrfull);
    remove(Kerrpolar);
}
return 0; }
}