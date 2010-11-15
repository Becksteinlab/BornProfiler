/* draw_membrane2.c                     09/02/08 *
 * draw_membrane2a.c                    11/13/10 * 
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
#include <assert.h>

typedef uint bool;

#define MAXLEN 100
#define TRUE  1
#define FALSE 0
#define NUMCOLS 3   /* number of columns in dx files, ABPS specific */

char *newname(const char *prefix, const char *infix, const char *suffix, const bool gzipped) {
  int l,m,n,p;
  char *s;
  static char default_suffix[] = ".dx";
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
  
  s = (char*)calloc(l+m+n+p+1, sizeof(char)); /* must be free in calling code! */
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
  char buf[MAXLEN]; 

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

int read_header(gzFile *in, int *dim_x, int *dim_y, int *dim_z,
		float *x0, float *y0, float *z0,
		float *dx, float *dy, float *dz,
		int *num_data) {
  char s[MAXLEN];
  float tmp;
  int tmp1;

  /* fixed format(!) APBS DX file header */
  gzgets(in,s,MAXLEN);
  gzgets(in,s,MAXLEN);
  gzgets(in,s,MAXLEN);
  gzgets(in,s,MAXLEN);   
  gzscanf(in, "%6s %1s %5s %13s %6s %i %i %i \n", s,s,s,s,s,dim_x,dim_y,dim_z);
  gzscanf(in, "%6s %f %f %f \n",s, x0, y0, z0);
  gzscanf(in, "%5s %f %f %f \n",s, dx, &tmp, &tmp);
  gzscanf(in, "%5s %f %f %f \n",s, &tmp, dy, &tmp);
  gzscanf(in, "%5s %f %f %f \n",s, &tmp, &tmp, dz);
  gzgets(in,s,MAXLEN);
  gzscanf(in, "%6s %i %5s %5s %4s %6s %4s %i %5s %i %4s %7s \n",s,&tmp1,s,s,s,s,s,&tmp1,s,num_data,s,s);
  return *num_data;
}

int read_data(gzFile *in, int num_data, float *x) {
  int n = 0;          /* data read from one line */
  int data_read = 0;  /* to check that we have everything */
  int idata = 0;      /* index of next data entry (sequential number) */
  int num_lines, num_last_entries;
  int i;

  /* should generate the read string using NUMCOLS 
     XXX: WILL NOT WORK unless NUMCOLS == 3 !!
  */
  assert(NUMCOLS == 3);

  num_last_entries = num_data % NUMCOLS; 
  num_lines = (num_data - num_last_entries)/NUMCOLS;  /* total data lines less one left in file */ 

  for (i=0; i<num_lines; i++)
    {
      n = gzscanf(in, "%f %f %f \n", &x[idata], &x[idata+1], &x[idata+2]);
      idata += NUMCOLS;
      data_read += n;
      if (n != NUMCOLS) {
	printf("ERROR [data line %d]: expected %d but got %d data\n", i+1, NUMCOLS, n);
      }
    }

  if (num_last_entries == 1) {
    n = gzscanf(in,"%f \n", &x[idata]);    /* reading in the last line */
    idata += 1;
  }
  else if (num_last_entries == 2) {
    n = gzscanf(in,"%f %f \n", &x[idata], &x[idata+1]);
    idata += 2;
  } 
  else {
    n = 0;
  }
  data_read += n;
  if (n != num_last_entries) {
    printf("ERROR [data line %d]: expected %d but got %d data\n", i+1, num_last_entries, n);
  }

  assert(idata == num_data); /* index now points to one beyond the end of the array */
  if (data_read != num_data) {
    printf("Expected %d data points but read %d: problem with file.\n", num_data, data_read);
    exit(EXIT_FAILURE);
  }
  return data_read;
}

void *xopen(const bool compression, char *filename, char *mode) {
  /* open filename in write mode either for gz or standard */
  void *stream; 
  if (compression)
    stream = gzopen(filename, mode);
  else
    stream = fopen(filename, mode);
  return stream;
}

int xclose(const bool compression, void *stream) {
  /* close filehandle for gz or standard */
  int status;
  if (compression)
    status = gzclose((gzFile*)stream);
  else
    status = fclose((FILE*)stream);
  return status;
}

/* note: it is not possible to write a "xprintf()" that passes a list
   of arguments down because there does not exists a gzprintf()
   function that accepts a va_list (and I am too lazy to write that
   one, too). See http://c-faq.com/varargs/handoff.html
*/

int write_header(const bool compression, void *out, char *name, float z_m0, float l_m,
		 int dim_x, int dim_y, int dim_z,
		 float x0, float y0, float z0,
		 float dx, float dy, float dz) {
  int status;
  const char fmt[] =   
    "# Data from draw_membrane2a.c\n"
    "#\n"
    "# %s with membrane: zmem = %4.2f, Lmem = %4.2f\n" //,name, z_m0, l_m
    "#\n"
    "object 1 class gridpositions counts %i %i %i \n" //, dim_x, dim_y, dim_z
    "origin %12.6E %12.6E %12.6E\n" //, x0, y0, z0
    "delta %12.6E %12.6E %12.6E\n" //, dx,0.000000E+00,0.000000E+00
    "delta %12.6E %12.6E %12.6E\n" //, 0.000000E+00,dy,0.000000E+00
    "delta %12.6E %12.6E %12.6E\n" //, 0.000000E+00,0.000000E+00,dz
    "object 2 class gridconnections counts %i %i %i\n" //, dim_x, dim_y, dim_z
    "object 3 class array type double rank 0 items %i data follows\n";  // num_data
  int num_data = dim_x * dim_y * dim_z;

  if (compression) 
    status = gzprintf((gzFile*)out, fmt, 
		      name, z_m0, l_m,
		      dim_x, dim_y, dim_z,
		      x0, y0, z0,
		      dx,    0E+00,  0E+00,
		      0E+00,    dy,  0E+00,
		      0E+00, 0E+00,     dz,
		      dim_x, dim_y, dim_z,
		      num_data);
  else
    status = fprintf((FILE*)out, fmt, 
		     name, z_m0, l_m,
		     dim_x, dim_y, dim_z,
		     x0, y0, z0,
		     dx,    0E+00,  0E+00,
		     0E+00,    dy,  0E+00,
		     0E+00, 0E+00,     dz,
		     dim_x, dim_y, dim_z,
		     num_data);
  return status;
}

#define WRITE_DATA(PRINTF, TYPE) do {		   \
   for (i=0; i<num_data; i++)                      \
     {	                                           \
       PRINTF((TYPE*)out, "%12.6E ", x[i]);        \
       if ((i+1) % NUMCOLS == 0)	           \
	 PRINTF((TYPE*)out, "\n"); /* break after NUMCOLS columns */ \
     }                                             \
   if (num_data % NUMCOLS != 0)                    \
     PRINTF((TYPE*)out, "\n");     /* finish incomplete last line */ \
  } while(0)

int write_data(const bool compression, void *out, int num_data, void *data) {
  float *x = (float *)data;     /* allows us to take int and write as float */
  int i;

  if (compression)
    WRITE_DATA(gzprintf, gzFile);
  else
    WRITE_DATA(fprintf, FILE);
  return i;
}


int write_attr_positions(const bool compression, void *stream) {
  int status;
  const char attr_positions[] = 
    "\n"
    "attribute \"dep\" string \"positions\"\n"
    "object \"regular positions regular connections\" class field\n"
    "component \"positions\" value 1\n"
    "component \"connections\" value 2\n"
    "component \"data\" value 3\n";
  if (compression)
    status = gzprintf((gzFile*)stream, attr_positions);
  else
    status = fprintf((FILE*)stream, attr_positions);
  return status;
}



/********************************************************************/
/* INPUT LOOKS LIKE:                                                */
/*    draw_membrane2a  infix  z_m0 l_m pdie V I R_m1 R_m0  gz       */
/*    draw_membrane2a  infix  z_m0 l_m pdie V I R_m1 R_m0           */
/********************************************************************/
void printhelp()
{
printf(         "* draw_membrane2a.c                    11/13/10 *\n"
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
int *map;
gzFile *out, *in; 
float *d_x, *d_y, *d_z;      
float *x_x, *x_y, *x_z;
float *y_x, *y_y, *y_z;
float *z_x, *z_y, *z_z;
float *kk, *cc;
float *x, *y, *z;
float tmp_x,tmp_y,tmp_z,dx,dy,dz,l_c_x,l_c_y,l_c_z,l_m;
float x0_p, y0_p, z0_p;
float x0_x, y0_x, z0_x;
float x0_y, y0_y, z0_y; 
float x0_z, y0_z, z0_z; 
float x0, y0, z0;
float V, I, sdie;
float z_m0, z_m1, R_m0, R_m1, R_x, R_y, R_z, R, pdie, mdie;
float R_temp;
char infix[MAXLEN];
char *file_name_x, *file_name_y, *file_name_z;
char *file_name_k, *file_name_c;
char *f1, *f2, *f3, *f4, *f5, *f6;
char ext[5]="m.dx";
bool compression = FALSE;

printf("----------------------------------------------------------------\n");
printf("draw_membrane2a -- (c) 2008 Michael Grabe [09/02/08]\n");
printf("                   (c) 2010 Oliver Beckstein, minor modifications [11/13/10]\n");
printf("Based on http://www.poissonboltzmann.org/apbs/examples/potentials-of-mean-force/the-polar-solvation-potential-of-mean-force-for-a-helix-in-a-dielectric-slab-membrane/draw_membrane2.c\n");
printf("----------------------------------------------------------------\n");

if (argc < 9) {
	printhelp();
	return 1;
}

strcpy(infix,argv[1]);
printf("Using hard-coded names with your infix to find files: infix=%s\n", infix);

if (argc == 10 && 0 == strcmp(argv[9], "gz")) {
  compression = TRUE;
  printf("Reading gzip-compressed dx files.\n");
} else {
  compression = FALSE;
  printf("Reading un-compressed dx files (default).\n");
}  

printf(">>> draw_membrane2a  %s  %s %s %s %s %s %s %s  %s\n",
       argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8],
       compression ? "gz" : "none");

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
read_header(in, &dim_x, &dim_y, &dim_z, &x0_x, &y0_x, &z0_x, &dx, &dy, &dz, &dim3);

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

read_data(in, dim3, d_x);

gzclose(in);
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
read_header(in, &dim_x, &dim_y, &dim_z, &x0_y, &y0_y, &z0_y, &dx, &dy, &dz, &dim3);

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

 read_data(in, dim3, d_y);

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
read_header(in, &dim_x, &dim_y, &dim_z, &x0_z, &y0_z, &z0_z, &dx, &dy, &dz, &dim3);

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
read_data(in, dim3, d_z);

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
read_header(in, &dim_x, &dim_y, &dim_z, &x0, &y0, &z0, &dx, &dy, &dz, &dim3);

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
 read_data(in, dim3, kk);

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

read_header(in, &dim_x, &dim_y, &dim_z, &x0, &y0, &z0, &dx, &dy, &dz, &dim3);

/* Read in the rest of the charge data */
read_data(in, dim3, cc);

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
out = xopen(compression, f1,"w");

/* MAKE THE X HEADER FILE */
 write_header(compression, out, "X-SHIFTED DIELECTRIC MAP", z_m0, l_m, 
	      dim_x, dim_y, dim_z, x0_x, y0_x, z0_x, dx ,dy, dz);

/* ADD THE X DATA */
 write_data(compression, out, dim3, d_x);

 write_attr_positions(compression, out);
 xclose(compression, out);
printf("Wrote %s.\n", f1);

/********************Y-DATA******************************/

/* give the file an "m" extension */
f2 = newname("diely", infix, ext, compression);
out = xopen(compression, f2,"w");
     
/* MAKE THE Y HEADER FILE */
write_header(compression, out, "Y-SHIFTED DIELECTRIC MAP", z_m0, l_m, 
	     dim_x, dim_y, dim_z, x0_y, y0_y, z0_y, dx ,dy, dz);

/* ADD THE Y DATA */
write_data(compression, out, dim3, d_y);

write_attr_positions(compression, out);
xclose(compression, out);
printf("Wrote %s.\n", f2);

/**********************Z-DATA*****************************/


/* give the file an "m" extension */
f3 = newname("dielz", infix, ext, compression);
out = xopen(compression, f3,"w");


/* MAKE THE Z HEADER FILE */
write_header(compression, out, "Z-SHIFTED DIELECTRIC MAP", z_m0, l_m, 
	     dim_x, dim_y, dim_z, x0_z, y0_z, z0_z, dx ,dy, dz);

/* ADD THE Z DATA */
write_data(compression, out, dim3, d_z);

write_attr_positions(compression, out);
xclose(compression, out);
printf("Wrote %s.\n", f3);

/*********************KAPPA******************************/

/* give the file an "m" extension */
f4 = newname("kappa", infix, ext, compression);
out = xopen(compression, f4,"w");

/* MAKE THE KAPPA HEADER FILE */
write_header(compression, out, "KAPPA MAP", z_m0, l_m, dim_x, dim_y, dim_z, x0, y0, z0, dx ,dy, dz);

/* ADD THE KAPPA DATA */
write_data(compression, out, dim3, kk);

write_attr_positions(compression, out);
xclose(compression, out);
printf("Wrote %s.\n", f4);

/********************CHARGE*******************************/

/* give the file an "m" extension */
f5 = newname("charge", infix, ext, compression);
out = xopen(compression, f5,"w");

/* MAKE THE CHARGE HEADER FILE */
write_header(compression, out, "CHARGE MAP", z_m0, l_m, dim_x, dim_y, dim_z, x0, y0, z0, dx ,dy, dz);

/* ADD THE CHARGE DATA */ 
write_data(compression, out, dim3, cc);

write_attr_positions(compression, out);
xclose(compression, out);
printf("Wrote %s.\n", f5);

/********************CHARGE CHANGE MAP*************************/

f6 = newname("change_map", infix, ext, compression);
out = xopen(compression, f6,"w");

/* MAKE THE CHARGE HEADER FILE */

write_header(compression, out, "CHARGE CHANGE MAP", z_m0, l_m, dim_x, dim_y, dim_z, x0, y0, z0, dx ,dy, dz);

/* ADD THE CHARGE CHANGE DATA */
 write_data(compression, out, dim3, map);  // int map is converted to float

write_attr_positions(compression, out);
xclose(compression, out);
printf("Wrote %s.\n", f6);

/***********************************************************/
free(x_x); free(y_x); free(z_x);
free(x_y); free(y_y); free(z_y); 
free(x_z); free(y_z); free(z_z); 
free(d_x); free(d_y); free(d_z);
free(x);   free(y);   free(z);
free(kk);  free(cc);  free(map);

free(f1); free(f2); free(f3); free(f4); free(f5); free(f6); /* calloc'ed by newname() :-p */

printf("Your files have been written.\n");
return 0;
}
