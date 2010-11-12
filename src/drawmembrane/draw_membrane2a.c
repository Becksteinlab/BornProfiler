/* draw_membrane2.c                     09/02/08 *
 * draw_membrane2a.c                    11/12/10 * 
 *-----------------------------------------------* 
 * By Michael Grabe                              *
 * This program takes the dielectric, kappa, and * 
 * charge  maps from APBS and accounts for the   *
 * membrane. the thickness and the bottom of the *
 * membrane must be given. we assume that the    *
 * membrane normal is along the z-axis. this     *
 * requires lining the protein along the z-axis  *
 * before running this program.                  *
 * We also output a charge_change_map that tells *
 * me which positions in the charge matrix were  *
 * edited by the addition of the membrane        *
 * NOTE: this program was changes in 2005 to     *
 * allow for conical channel geometries. An      *
 * extra command line argument was added that    *
 * will cause a fault if run without it from     *
 * older scripts pre 2005.                       *
 *                                               * 
 * INPUTS:                                       *
 * thses all come at the command line:           *
 * infix  - infix is used to construct map names *
 * z_m0   - bottom of the membrane               *
 * l_m    - length of the membrane               *
 * pdie   - protein dielectric constant          *
 * V      - cytoplasmic potential (kT/e)         *
 * I      - molar conc. of one salt-species      *
 * R_m1   - excl. radius at top of  membrane     *
 * R_m0   - excl. radius at  bottom of membrane  *
 *                                               * 
 * OUTPUTS:                                      *
 *   maps                                        *
 *-----------------------------------------------*

 2010-11-11  Oliver Beckstein
             changed naming (provide the infix) and added more diagnostics
 2010-11-12  Oliver Beckstein
             added reading/writing of gz-compressed files
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <zlib.h>

typedef uint bool;

#define MAXLEN 100
#define TRUE  1
#define FALSE 0

char *newname(char *prefix, char *infix, char *suffix, bool gzipped) {
  int l,m,n,p;
  char *s;
  static char default_suffix[4] = ".dx";
  char compression_suffix[4];  /* ".gz" or "" */

  if (NULL == suffix) {
    suffix = default_suffix;
  }

  if (gzipped) {
    strcpy(compression_suffix, ".gz");
  } else {
    strcpy(compression_suffix, "");
  }

  l = strlen(prefix);
  m = strlen(infix);
  n = strlen(suffix);
  p = strlen(compression_suffix);
  
  s = (char*)calloc(l+m+n+p+1, sizeof(char)); /* not freed ever... */
  if (NULL == s) {
    printf("newname: Failed to allocate string.");
    exit(EXIT_FAILURE);
  }
  strncpy(s, prefix, l);
  strncat(s, infix, m);
  strncat(s, suffix, n);
  strncat(s, compression_suffix, n);
  return s;
}


int gzscanf(gzFile *stream, const char *fmt, ...) {
  /* read one line from stream (up to newline) and parse with sscanf */
  va_list args;
  va_start(args, fmt);
  int n;
  static char buf[MAXLEN]; 

  if (NULL == gzgets(stream, buf, MAXLEN)) {
    printf("gzscanf: Failed to read line from gz file.\n");
    exit(EXIT_FAILURE);
  }
  n = vsscanf(buf, fmt, args);
  va_end(args);
  if (0 == n) {
    printf("gzscanf: Failed to parse line: %s\n", buf);
    exit(EXIT_FAILURE);
  }
  return n;
}

/********************************************************************/
/* INPUT LOOKS LIKE:                                                */
/*    draw_membrane2a  infix  z_m0 l_m pdie V I R_m1 R_m0  gz       */
/*    draw_membrane2a  infix  z_m0 l_m pdie V I R_m1 R_m0           */
/********************************************************************/
void printhelp()
{
printf(         "* draw_membrane2a.c                    11/12/10 *\n"
		"*-----------------------------------------------*\n"
		"* By Michael Grabe                              *\n"
                "* (minor modifications by Oliver Beckstein)     *\n"
		"* This program takes the gz-compressed          *\n"
		"* dielectric, kappa, and                        *\n" 
		"* charge  maps from APBS and accounts for the   *\n"
		"* membrane. the thickness and the bottom of the *\n"
		"* membrane must be given. we assume that the    *\n"
		"* membrane normal is along the z-axis. this     *\n"
		"* requires lining the protein along the z-axis  *\n"
		"* before running this program.                  *\n"
		"* We also output a charge_change_map that tells *\n"
		"* me which positions in the charge matrix were  *\n"
		"* edited by the addition of the membrane        *\n"
		"* NOTE: this program was changes in 2005 to     *\n"
		"* allow for conical channel geometries. An      *\n"
		"* extra command line argument was added that    *\n"
		"* will cause a fault if run without it from     *\n"
		"* older scripts pre 2005.                       *\n"
		"*                                               *\n" 
		"* INPUTS:                                       *\n"
		"* these all come at the command line:           *\n"
		"* infix  - all maps are constructed as          *\n"
                "*          <name><infix>.dx[.gz] where <name> is*\n"
                "*          hard-coded (dielx,diely,dielz,kappa  *\n"
                "*          charge); see also OUTPUTS.           *\n"
		"* z_m0   - bottom of the membrane               *\n"
		"* l_m    - length of the membrane               *\n"
		"* pdie   - protein dielectric constant          *\n"
		"* V      - cytoplasmic potential (kT/e)         *\n"
		"* I      - molar conc. of one salt-species      *\n"
		"* R_m1   - excl. radius at top of  membrane     *\n"
		"* R_m0   - excl. radius at  bottom of membrane  *\n"
                "* gz     - set to 'gz' if files are compressed  *\n"
                "*          leave out or 'none' otherwise        *\n"
		"*                                               *\n" 
		"* OUTPUTS:                                      *\n"
		"*   maps - names are <name><infix>m.dx[.gz]     *\n"
		"*-----------------------------------------------*\n\n"
		"********************************************************************\n"
		"* INPUT LOOKS LIKE:                                                *\n"
		"* ./draw_membrane2a infix  z_m0 l_m pdie V I R_m1 R_m0  [gz]       *\n"
		"********************************************************************\n\n"
	  );
}


int main(int argc, char *argv[])
{

int dim_x,dim_y,dim_z,dim3,i,j,k,cnt;
int tmp1,tmp2,tmp3,tmp4;
int l, *map;
gzFile *out, *in; 
float *d_x, *d_y, *d_z;      
float *x_x, *x_y, *x_z;
float *y_x, *y_y, *y_z;
float *z_x, *z_y, *z_z;
float *kk, *cc;
float *x, *y, *z;
float tmp,tmp_x,tmp_y,tmp_z,dx,dy,dz,l_c_x,l_c_y,l_c_z,l_m;
float x0_p, y0_p, z0_p;
float x0_x, y0_x, z0_x;
float x0_y, y0_y, z0_y; 
float x0_z, y0_z, z0_z; 
float x0, y0, z0;
float V, I, sdie;
float z_m0, z_m1, R_m0, R_m1, R_x, R_y, R_z, R, pdie, mdie;
float R_temp;
char infix[MAXLEN];
char s[MAXLEN];
char *file_name_x, *file_name_y, *file_name_z;
char *file_name_k, *file_name_c;
char *f1, *f2, *f3, *f4, *f5, *f6;
char ext[5]="m.dx";
bool compression = FALSE;

if (argc < 9) {
	printhelp();
	return 1;
}

printf("----------------------------------------------------------------\n");
printf("draw_membrane2a -- (c) 2008 Michael Grabe [09/02/08]\n");
printf("                   (c) 2010 Oliver Beckstein, minor modifications [11/12/10]\n");
printf("Based on http://www.poissonboltzmann.org/apbs/examples/potentials-of-mean-force/the-polar-solvation-potential-of-mean-force-for-a-helix-in-a-dielectric-slab-membrane/draw_membrane2.c\n");
printf("----------------------------------------------------------------\n");

strcpy(infix,argv[1]);
printf("Using hard-coded names with your infix to find files: infix=%s\n", infix);

if (argc == 10 && 0 == strcmp(argv[9], "gz")) {
  compression = TRUE;
  printf("Reading gzip-compressed dx files.\n");
} else {
  compression = FALSE;
  printf("Reading un-compressed dx files (default).\n");
}  

/* Find the x-shifted dielectric map 
   Construct the name as <basename><infix><suffix>

   suffix == NULL --> use default = ".dx"
 */
file_name_x = newname("dielx", infix, NULL, compression);

/* Find the y-shifted dielectric map */
file_name_y = newname("diely", infix, NULL, compression);

/* Find the z-shifted dielectric map */
file_name_z = newname("dielz", infix, NULL, compression);

/* Find the kappa map */
file_name_k = newname("kappa", infix, NULL, compression);

/* Find the charge map */
file_name_c = newname("charge", infix, NULL, compression);

z_m0=atof(argv[2]);
l_m=atof(argv[3]);
pdie=atof(argv[4]);
V=atof(argv[5]);
I=atof(argv[6]);

R_m1=atof(argv[7]);
R_m0=atof(argv[8]);

z_m1=z_m0+l_m;   /* top of the membrane */
mdie = 2.0;    /* watch out for this it used to be 10.0 */ 
sdie = 80.0;
/*****************************************************/
/* read in the x-shifted dielectric data             */
/*****************************************************/

in = gzopen(file_name_x,"r");
if (in == NULL) {
	printhelp();
	printf("Make sure %s exists in current directory!!!\n\n", file_name_x);
	return 1;
}
printf("Reading %s...\n", file_name_x);

/* First read the header */

gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);   
gzscanf(in, "%6s %1s %5s %13s %6s %i %i %i \n", s,s,s,s,s,&dim_x,&dim_y,&dim_z);
gzscanf(in, "%6s %f %f %f \n",s, &x0_x, &y0_x, &z0_x);
gzscanf(in, "%5s %f %f %f \n",s, &dx, &tmp, &tmp);
gzscanf(in, "%5s %f %f %f \n",s, &tmp, &dy, &tmp);
gzscanf(in, "%5s %f %f %f \n",s, &tmp, &tmp, &dz);
gzgets(in,s,MAXLEN);
gzscanf(in, "%6s %i %5s %5s %4s %6s %4s %i %5s %i %4s %7s \n",s,&tmp1,s,s,s,s,s,&tmp1,s,&dim3,s,s);


/* assign the memory to the arrays */

x_x= (float *) calloc(dim_x+1,sizeof(float));
y_x= (float *) calloc(dim_y+1,sizeof(float));
z_x= (float *) calloc(dim_z+1,sizeof(float));
x_y= (float *) calloc(dim_x+1,sizeof(float));
y_y= (float *) calloc(dim_y+1,sizeof(float));
z_y= (float *) calloc(dim_z+1,sizeof(float));
x_z= (float *) calloc(dim_x+1,sizeof(float));
y_z= (float *) calloc(dim_y+1,sizeof(float));
z_z= (float *) calloc(dim_z+1,sizeof(float));
d_x= (float *) calloc(dim_x*dim_y*dim_z+1,sizeof(float));
d_y= (float *) calloc(dim_x*dim_y*dim_z+1,sizeof(float));
d_z= (float *) calloc(dim_x*dim_y*dim_z+1,sizeof(float));
/* Now the Kappa and charge Arrays */
x= (float *) calloc(dim_x+1,sizeof(float));
y= (float *) calloc(dim_y+1,sizeof(float));
z= (float *) calloc(dim_z+1,sizeof(float));
kk= (float *) calloc(dim_x*dim_y*dim_z+1,sizeof(float));
cc= (float *) calloc(dim_x*dim_y*dim_z+1,sizeof(float));
map= (int *) calloc(dim_x*dim_y*dim_z+1,sizeof(int));

/* initialize x,y,z, and diel vectors */

l_c_x=dim_x*dx;
l_c_y=dim_y*dx;
l_c_z=dim_z*dx;

tmp_x=x0_x;
tmp_y=y0_x;
tmp_z=z0_x;

for (i=1; i <= dim_x; ++i)
{
x_x[i]=tmp_x;
tmp_x+=dx;
}

for (i=1; i <= dim_y; ++i)
{
y_x[i]=tmp_y;
tmp_y+=dy;
}

for (i=1; i <= dim_z; ++i)
{
z_x[i]=tmp_z;
tmp_z+=dz;
}

/* Read in the rest of the dielectric data */

tmp1 = fmod(dim3, 3); 
tmp2 =(dim3-tmp1)/3;    /* total lines less one left in file */ 
tmp3 =1;

for (i=1; i <= tmp2; ++i)
{
gzscanf(in,"%f %f %f \n", &d_x[tmp3], &d_x[tmp3+1], &d_x[tmp3+2]); 
tmp3+=3;
}

if (tmp1 == 1)
gzscanf(in,"%f \n", &d_x[tmp3]);    /* reading in the last line */
else if (tmp1 == 2)
gzscanf(in,"%f %f \n", &d_x[tmp3], &d_x[tmp3+1]);

/*****************************************************/

/****************************************************/
/* Construct the protein cooridinate center based   */
/* on the x-dielectric map.                         */
/* This will be used to determine where to add      */
/* membrane.                                        */ 
/****************************************************/

x0_p=x0_x+l_c_x/2-dx/2;  /* this is the shift term that */ 
                         /* moves half-step off grid    */ 
y0_p=y0_x+l_c_y/2;
z0_p=z0_x+l_c_z/2;

gzclose(in);

/*****************************************************/
/* read in the y-shifted dielectric data             */
/*****************************************************/

in = gzopen(file_name_y,"r");
if (in == NULL) {
   printf("File name %s not found.\n", file_name_y);
   return 1;
}
printf("Reading %s...\n", file_name_y);

/* First read the header */

gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzscanf(in, "%6s %f %f %f \n",s, &x0_y, &y0_y, &z0_y);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzscanf(in, "%6s %i %5s %5s %4s %6s %4s %i %5s %i %4s %7s \n",s,&tmp1,s,s,s,s,s,&tmp1,s, &tmp1,s,s); 
 
/* initialize x,y,z, and diel vectors */

tmp_x=x0_y;
tmp_y=y0_y;
tmp_z=z0_y;

for (i=1; i <= dim_x; ++i) {
x_y[i]=tmp_x;
tmp_x+=dx;
}

for (i=1; i <= dim_y; ++i)
{
y_y[i]=tmp_y;
tmp_y+=dy;
}

for (i=1; i <= dim_z; ++i)
{
z_y[i]=tmp_z;
tmp_z+=dz;
}

/* Read in the rest of the dielectric data */

tmp3 =1;

for (i=1; i <= tmp2; ++i)
{
gzscanf(in,"%f %f %f \n", &d_y[tmp3], &d_y[tmp3+1], &d_y[tmp3+2]);
tmp3+=3;
}

if (tmp1 == 1)
gzscanf(in,"%f \n", &d_y[tmp3]);    /* reading in the last line */
else if (tmp1 == 2)
gzscanf(in,"%f %f \n", &d_y[tmp3], &d_y[tmp3+1]);

gzclose(in);

/*****************************************************/
/* read in the z-shifted dielectric data             */
/*****************************************************/

in = gzopen(file_name_z,"r");
if (in == NULL) {
   printf("File name %s not found.\n", file_name_z);
   return 1;
}
printf("Reading %s...\n", file_name_z);

/* First read the header */

gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzscanf(in, "%6s %f %f %f \n",s, &x0_z, &y0_z, &z0_z);
gzgets(in,s,MAXLEN);                                                 
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzscanf(in, "%6s %i %5s %5s %4s %6s %4s %i %5s %i %4s %7s \n",s,&tmp1,s,s,s,s,s,&tmp1,s, &tmp1,s,s);

/* initialize x,y,z, and diel vectors */

tmp_x=x0_z;
tmp_y=y0_z;
tmp_z=z0_z;


for (i=1; i <= dim_x; ++i)
{
x_z[i]=tmp_x;
tmp_x+=dx;
}

for (i=1; i <= dim_y; ++i)
{
y_z[i]=tmp_y;
tmp_y+=dy;
}

for (i=1; i <= dim_z; ++i)
{
z_z[i]=tmp_z;
tmp_z+=dz;
}


/* Read in the rest of the dielectric data */

tmp3 =1;

for (i=1; i <= tmp2; ++i)
{ 
gzscanf(in,"%f %f %f \n", &d_z[tmp3], &d_z[tmp3+1], &d_z[tmp3+2]);
tmp3+=3;
}

if (tmp1 == 1) 
gzscanf(in,"%f \n", &d_z[tmp3]);    /* reading in the last line */
else if (tmp1 == 2)
gzscanf(in,"%f %f \n", &d_z[tmp3], &d_z[tmp3+1]);

gzclose(in);

/*****************************************************/

/*****************************************************/
/* read in the kappa data                            */
/*****************************************************/

in = gzopen(file_name_k,"r");
if (in == NULL) {
   printf("File name %s not found.\n", file_name_k);
   return 1;
}
printf("Reading %s...\n", file_name_k);

/* First read the header */

gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzscanf(in, "%6s %f %f %f \n",s, &x0, &y0, &z0);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);

/* initialize x,y,z, and kappa vectors */

tmp_x=x0;
tmp_y=y0;
tmp_z=z0;


for (i=1; i <= dim_x; ++i)
{
x[i]=tmp_x;
tmp_x+=dx;
}

for (i=1; i <= dim_y; ++i)
{
y[i]=tmp_y;
tmp_y+=dy;
}

for (i=1; i <= dim_z; ++i)
{
z[i]=tmp_z;
tmp_z+=dz;
}

/* Read in the rest of the Kappa data */


tmp3 =1;

for (i=1; i <= tmp2; ++i)
{
gzscanf(in,"%f %f %f \n", &kk[tmp3], &kk[tmp3+1], &kk[tmp3+2]);
tmp3+=3;
}

if (tmp1 == 1)
gzscanf(in,"%f \n", &kk[tmp3]);    /* reading in the last line */
else if (tmp1 == 2)
gzscanf(in,"%f %f \n", &kk[tmp3], &kk[tmp3+1]);

gzclose(in);

/*****************************************************/

/*****************************************************/
/* read in the charge data                           */
/*****************************************************/

in = gzopen(file_name_c,"r");
if (in == NULL) {
   printf("File name %s not found.\n", file_name_c);
   return 1;
}
printf("Reading %s...\n", file_name_c);

/* First read the header */

gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);
gzgets(in,s,MAXLEN);

/* Read in the rest of the charge data */

tmp3 =1;

for (i=1; i <= tmp2; ++i)
{
gzscanf(in,"%f %f %f \n", &cc[tmp3], &cc[tmp3+1], &cc[tmp3+2]);
tmp3+=3;
}

if (tmp1 == 1)
gzscanf(in,"%f \n", &cc[tmp3]);    /* reading in the last line */
else if (tmp1 == 2)
gzscanf(in,"%f %f \n", &cc[tmp3], &cc[tmp3+1]);

gzclose(in);

/*****************************************************/
/* MANIPULATE THE DATA BY ADDING THE MEMBRANE        */      
/*****************************************************/


/******************************************************/
/* set up the vector                                  */
/******************************************************/


cnt=1;

for (k=1; k <= dim_x; ++k)  /* loop over z */
{
        for (j=1; j <= dim_y; ++j)  /* loop over y */
	{	
		for (i=1; i <= dim_z; ++i)  /* loop over x */
	        {	
                        R_x = sqrt((x_x[k]-x0_p)*(x_x[k]-x0_p) + (y_x[j]-y0_p)*(y_x[j]-y0_p));	
		        R_temp = (R_m1*(z_x[i]-z_m0) - R_m0*(z_x[i]-z_m1))/(z_m1 - z_m0);  	

                        if (z_x[i] <= z_m1 && z_x[i] >= z_m0 && R_x > R_temp && d_x[cnt] > pdie+0.05) 
		        {	
                        d_x[cnt] = mdie;   /* bilayer dielectric constant */
                        }

                        R_y = sqrt((x_y[k]-x0_p)*(x_y[k]-x0_p) + (y_y[j]-y0_p)*(y_y[j]-y0_p));
                        R_temp = (R_m1*(z_y[i]-z_m0) - R_m0*(z_y[i]-z_m1))/(z_m1 - z_m0);
 
                        if (z_y[i] <= z_m1 && z_y[i] >= z_m0 && R_y > R_temp && d_y[cnt] > pdie+0.05)
                        {      
                        d_y[cnt] = mdie;   /* bilayer dielectric constant */
                        }

                        R_z = sqrt((x_z[k]-x0_p)*(x_z[k]-x0_p) + (y_z[j]-y0_p)*(y_z[j]-y0_p));
                        R_temp = (R_m1*(z_z[i]-z_m0) - R_m0*(z_z[i]-z_m1))/(z_m1 - z_m0);

                        if (z_z[i] <= z_m1 && z_z[i] >= z_m0 && R_z > R_temp && d_z[cnt] > pdie+0.05)
                        {      
                        d_z[cnt] = mdie;   /* bilayer dielectric constant */
                        }

                        R = sqrt((x[k]-x0_p)*(x[k]-x0_p) + (y[j]-y0_p)*(y[j]-y0_p));

                        if (z[i] <= z_m0 && kk[cnt] != 0.0)
                        {
                        /* charge for mem V */
                        /* see my notes for this expression */
                        cc[cnt] = 0.0012045*I*V;
                        /* update the change map */
                        map[cnt] = 1;
                        }
                        else
                        {
                        /* position was not changed */
                        map[cnt] = 0;
                        }
                      
                        R_temp = (R_m1*(z[i]-z_m0) - R_m0*(z[i]-z_m1))/(z_m1 - z_m0);

                        if (z[i] <= z_m1 && z[i] >= z_m0 && R > R_temp)
                        {
                        kk[cnt] = 0.0;   /* Zero ion accessibility */
                        }

                        ++cnt;

  	        }			
	}

}

/********************************************************/
/* now we must save diel as a text document             */
/* in the proper 3 column output                        */
/********************************************************/



/********************************************************/

/* add the "m" extension to the file */
 f1 = newname("dielx", infix, ext, compression);
out = gzopen(f1,"w");

/* MAKE THE X HEADER FILE */

gzprintf(out, "# Data from draw_membrane.c \n");
gzprintf(out, "# \n");
gzprintf(out, "# X-SHIFTED DIELECTRIC MAP with membrane: zmem = %4.2f, Lmem = %4.2f \n",z_m0, l_m);
gzprintf(out, "# \n");
gzprintf(out, "object 1 class gridpositions counts %i %i %i \n", dim_x, dim_y, dim_z);
gzprintf(out, "origin %12.6E %12.6E %12.6E \n", x0_x, y0_x, z0_x);
gzprintf(out, "delta %12.6E %12.6E %12.6E \n", dx,0.000000E+00,0.000000E+00);
gzprintf(out, "delta %12.6E %12.6E %12.6E \n", 0.000000E+00,dy,0.000000E+00);
gzprintf(out, "delta %12.6E %12.6E %12.6E \n", 0.000000E+00,0.000000E+00,dz);
gzprintf(out, "object 2 class gridconnections counts %i %i %i\n", dim_x, dim_y, dim_z);
gzprintf(out, "object 3 class array type double rank 0 items %i data follows\n", dim_x*dim_y*dim_z);

/* ADD THE X DATA */

cnt=1;

for (i=1; i <= tmp2  ; ++i)
{	
gzprintf(out, "%12.6E %12.6E %12.6E \n", d_x[cnt], d_x[cnt+1], d_x[cnt+2]);
cnt=cnt+3;
}


if (tmp1 == 1)
gzprintf(out,"%12.6E \n", d_x[cnt]);    /* saving in the last line */
else if (tmp1 == 2)
gzprintf(out,"%12.6E %12.6E \n", d_x[cnt], d_x[cnt+1]);

gzclose(out);
printf("Wrote %s.\n", f1);

/********************Y-DATA******************************/

/* give the file an "m" extension */
f2 = newname("diely", infix, ext, compression);
out = gzopen(f2,"w");
     
/* MAKE THE Y HEADER FILE */

gzprintf(out, "# Data from draw_membrane.c \n");
gzprintf(out, "# \n");
gzprintf(out, "# Y-SHIFTED DIELECTRIC MAP with membrane: zmem = %4.2f, Lmem = %4.2f \n",z_m0, l_m);
gzprintf(out, "# \n");
gzprintf(out, "object 1 class gridpositions counts %i %i %i \n", dim_x, dim_y, dim_z);
gzprintf(out, "origin %12.6E %12.6E %12.6E \n", x0_y, y0_y, z0_y);
gzprintf(out, "delta %12.6E %12.6E %12.6E \n", dx,0.000000E+00,0.000000E+00);
gzprintf(out, "delta %12.6E %12.6E %12.6E \n", 0.000000E+00,dy,0.000000E+00);
gzprintf(out, "delta %12.6E %12.6E %12.6E \n", 0.000000E+00,0.000000E+00,dz);
gzprintf(out, "object 2 class gridconnections counts %i %i %i\n", dim_x, dim_y, dim_z)
;
gzprintf(out, "object 3 class array type double rank 0 items %i data follows\n", dim_x*dim_y*dim_z);

/* ADD THE Y DATA */

cnt=1;

for (i=1; i <= tmp2  ; ++i)
{      
gzprintf(out, "%12.6E %12.6E %12.6E \n", d_y[cnt], d_y[cnt+1], d_y[cnt+2]);
cnt=cnt+3;
}

if (tmp1 == 1)
gzprintf(out,"%12.6E \n", d_y[cnt]);    /* saving in the last line */ 
else if (tmp1 == 2)
gzprintf(out,"%12.6E %12.6E \n", d_y[cnt], d_y[cnt+1]);

gzclose(out);
printf("Wrote %s.\n", f2);

/**********************Z-DATA*****************************/


/* give the file an "m" extension */
f3 = newname("dielz", infix, ext, compression);
out = gzopen(f3,"w");


/* MAKE THE Z HEADER FILE */

gzprintf(out, "# Data from draw_membrane.c \n");
gzprintf(out, "# \n");
gzprintf(out, "# Z-SHIFTED DIELECTRIC MAP with membrane: zmem = %4.2f, Lmem = %4.2f \n",z_m0, l_m);
gzprintf(out, "# \n");
gzprintf(out, "object 1 class gridpositions counts %i %i %i \n", dim_x, dim_y, dim_z);
gzprintf(out, "origin %12.6E %12.6E %12.6E \n", x0_z, y0_z, z0_z);
gzprintf(out, "delta %12.6E %12.6E %12.6E \n", dx,0.000000E+00,0.000000E+00);
gzprintf(out, "delta %12.6E %12.6E %12.6E \n", 0.000000E+00,dy,0.000000E+00);
gzprintf(out, "delta %12.6E %12.6E %12.6E \n", 0.000000E+00,0.000000E+00,dz);
gzprintf(out, "object 2 class gridconnections counts %i %i %i\n", dim_x, dim_y, dim_z);
gzprintf(out, "object 3 class array type double rank 0 items %i data follows\n", dim_x*dim_y*dim_z);

/* ADD THE Z DATA */

cnt=1;

for (i=1; i <= tmp2  ; ++i)
{      
gzprintf(out, "%12.6E %12.6E %12.6E \n", d_z[cnt], d_z[cnt+1], d_z[cnt+2]);
cnt=cnt+3;
}


if (tmp1 == 1)
gzprintf(out,"%12.6E \n", d_z[cnt]);    /* saving in the last line */ 
else if (tmp1 == 2)
gzprintf(out,"%12.6E %12.6E \n", d_z[cnt], d_z[cnt+1]);

gzclose(out);
printf("Wrote %s.\n", f3);

/*********************KAPPA******************************/

/* give the file an "m" extension */
f4 = newname("kappa", infix, ext, compression);
out = gzopen(f4,"w");


/* MAKE THE KAPPA HEADER FILE */

gzprintf(out, "# Data from draw_membrane.c \n");
gzprintf(out, "# \n");
gzprintf(out, "# KAPPA MAP with membrane: zmem = %4.2f, Lmem = %4.2f \n",z_m0, l_m);
gzprintf(out, "# \n");
gzprintf(out, "object 1 class gridpositions counts %i %i %i \n", dim_x, dim_y, dim_z);
gzprintf(out, "origin %12.6E %12.6E %12.6E \n", x0, y0, z0);
gzprintf(out, "delta %12.6E %12.6E %12.6E \n", dx,0.000000E+00,0.000000E+00);
gzprintf(out, "delta %12.6E %12.6E %12.6E \n", 0.000000E+00,dy,0.000000E+00);
gzprintf(out, "delta %12.6E %12.6E %12.6E \n", 0.000000E+00,0.000000E+00,dz);
gzprintf(out, "object 2 class gridconnections counts %i %i %i\n", dim_x, dim_y, dim_z);
gzprintf(out, "object 3 class array type double rank 0 items %i data follows\n", dim_x*dim_y*dim_z);

/* ADD THE KAPPA DATA */

cnt=1;

for (i=1; i <= tmp2  ; ++i)
{
gzprintf(out, "%12.6E %12.6E %12.6E \n", kk[cnt], kk[cnt+1], kk[cnt+2]);
cnt=cnt+3;
}

if (tmp1 == 1)
gzprintf(out,"%12.6E \n", kk[cnt]);    /* saving in the last line */
else if (tmp1 == 2)
gzprintf(out,"%12.6E %12.6E \n", kk[cnt], kk[cnt+1]);

gzclose(out);
printf("Wrote %s.\n", f4);

/********************CHARGE*******************************/

/* give the file an "m" extension */
f5 = newname("charge", infix, ext, compression);
out = gzopen(f5,"w");

/* MAKE THE CHARGE HEADER FILE */

gzprintf(out, "# Data from draw_membrane.c \n");
gzprintf(out, "# \n");
gzprintf(out, "# CHARGE MAP with membrane: zmem = %4.2f, Lmem = %4.2f \n",z_m0, l_m);
gzprintf(out, "# \n");
gzprintf(out, "object 1 class gridpositions counts %i %i %i \n", dim_x, dim_y, dim_z);
gzprintf(out, "origin %12.6E %12.6E %12.6E \n", x0, y0, z0);      
gzprintf(out, "delta %12.6E %12.6E %12.6E \n", dx,0.000000E+00,0.000000E+00);
gzprintf(out, "delta %12.6E %12.6E %12.6E \n", 0.000000E+00,dy,0.000000E+00); 
gzprintf(out, "delta %12.6E %12.6E %12.6E \n", 0.000000E+00,0.000000E+00,dz);
gzprintf(out, "object 2 class gridconnections counts %i %i %i\n", dim_x, dim_y, dim_z);
gzprintf(out, "object 3 class array type double rank 0 items %i data follows\n", dim_x*dim_y*dim_z);

/* ADD THE CHARGE DATA */ 

cnt=1;

for (i=1; i <= tmp2  ; ++i)
{
gzprintf(out, "%12.6E %12.6E %12.6E \n", cc[cnt], cc[cnt+1], cc[cnt+2]);   
cnt=cnt+3;
}

if (tmp1 == 1)
gzprintf(out,"%12.6E \n", cc[cnt]);    /* saving in the last line */ 
else if (tmp1 == 2)
gzprintf(out,"%12.6E %12.6E \n", cc[cnt], cc[cnt+1]);  

gzclose(out);
printf("Wrote %s.\n", f5);

/********************CHARGE CHANGE MAP*************************/

f6 = newname("change_map", infix, ext, compression);
out = gzopen(f6,"w");

/* MAKE THE CHARGE HEADER FILE */

gzprintf(out, "# Data from draw_membrane.c \n");
gzprintf(out, "# \n");
gzprintf(out, "# CHARGE CHANGE MAP with membrane: zmem = %4.2f, Lmem = %4.2f \n",z_m0, l_m);
gzprintf(out, "# \n");
gzprintf(out, "object 1 class gridpositions counts %i %i %i \n", dim_x, dim_y, dim_z);
gzprintf(out, "origin %12.6E %12.6E %12.6E \n", x0, y0, z0);
gzprintf(out, "delta %12.6E %12.6E %12.6E \n", dx,0.000000E+00,0.000000E+00);
gzprintf(out, "delta %12.6E %12.6E %12.6E \n", 0.000000E+00,dy,0.000000E+00);
gzprintf(out, "delta %12.6E %12.6E %12.6E \n", 0.000000E+00,0.000000E+00,dz);
gzprintf(out, "object 2 class gridconnections counts %i %i %i\n", dim_x, dim_y, dim_z);
gzprintf(out, "object 3 class array type double rank 0 items %i data follows\n",
dim_x*dim_y*dim_z);

/* ADD THE CHARGE CHANGE DATA */

cnt=1;

for (i=1; i <= tmp2  ; ++i)
{
gzprintf(out, "%i %i %i \n", map[cnt], map[cnt+1], map[cnt+2]);
cnt=cnt+3;
}

if (tmp1 == 1)
gzprintf(out,"%i \n", map[cnt]);    /* saving in the last line */
else if (tmp1 == 2)
gzprintf(out,"%i %i \n", map[cnt], map[cnt+1]);

gzclose(out);
printf("Wrote %s.\n", f6);

/***********************************************************/
free(x_x);
free(y_x);
free(z_x);
free(x_y);
free(y_y);
free(z_y); 
free(x_z);
free(y_z);
free(z_z); 
free(x);
free(y);
free(z);

printf("Your files have been written.\n");
return 0;
}
