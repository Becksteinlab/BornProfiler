/* draw_membrane2.c                     09/02/08 *
 * draw_membrane2a.c                    11/22/10 * 
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
 * OB's changes completely break input as we are now using standard option
 * processing. 

 2010-11-11  Oliver Beckstein
             changed naming (provide the infix) and added more diagnostics
 2010-11-12  Oliver Beckstein
             added reading/writing of gz-compressed files
 2010-11-22  OB: switched to opt processing (completely breaks interface
             but allows setting of sdie, mdie and pdie)
 2010-12-01  OB: added cdie cylinder dielectric
             headgroups (from draw_membrane4.c) and pretty-print geometry overview
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <zlib.h>
#include <assert.h>


#define SQR(x) ((x)*(x))

#define MAXLEN 100
#define TRUE  1
#define FALSE 0
#define NUMCOLS 3   /* number of columns in dx files, ABPS specific */

typedef uint bool;
typedef struct {
  float z_m0; /* bottom of membrane */
  float z_m1; /* top of membrane */
  float z_h0; /* upper boundary of lower headgroups, z_m0 + l_h */
  float z_h1; /* lower boundary of upper headgroups, z_m1 - l_h */
  float pdie; /* protein dielectric (must be SAME as in APBS!) */
  float mdie; /* membrane dielectric */
  float cdie; /* channel dielectric (defaults to solvent dielectric) */
  float idie; /* headgroup dielectric (defaults to mdie) */
  float R_m0; /* exclusion radii */
  float R_m1;
  float x0_p; /* x and y of protein */
  float y0_p;
  float z0_p;
} t_membrane;

char *newname(const char *prefix, const char *infix, const char *suffix, const bool gzipped);
int gzscanf(gzFile *stream, const char *fmt, ...);
void *xopen(const bool compression, char *filename, char *mode);
int xclose(const bool compression, void *stream);
int read_header(gzFile *in, int *dim_x, int *dim_y, int *dim_z,
		float *x0, float *y0, float *z0,
		float *dx, float *dy, float *dz,
		int *num_data);
int read_data(gzFile *in, int num_data, float *x);
int write_header(const bool compression, void *out, char *name, float z_m0, float l_m,
		 int dim_x, int dim_y, int dim_z,
		 float x0, float y0, float z0,
		 float dx, float dy, float dz);
int write_data(const bool compression, void *out, int num_data, void *data);
int write_attr_positions(const bool compression, void *stream);
float draw_diel(const t_membrane *M, float *d, const float x, const float y, const float z);
void print_help();
void print_membrane(const t_membrane *M);
void print_exclusionzone(const t_membrane *M);
int get_argv_int(char **argv, const int index);
int get_argv_float(char **argv, const int index);
char *get_argv_str(char **argv, const int index);

/**************************************************************/

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
  
  s = calloc(l+m+n+p+1, sizeof(char)); /* must be free'ed in calling code! */
  if (NULL == s) {
    fprintf(stderr, "newname: Failed to allocate string.");
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
    fprintf(stderr, "gzscanf: Failed to read line from gz file.\n");
    exit(EXIT_FAILURE);
  }
  n = vsscanf(buf, fmt, args);
  va_end(args);
  if (0 == n) {
    fprintf(stderr, "gzscanf: Failed to parse line: %s\n", buf);
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
	fprintf(stderr, "ERROR [data line %d]: expected %d but got %d data\n", i+1, NUMCOLS, n);
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
    fprintf(stderr, "ERROR [data line %d]: expected %d but got %d data\n", i+1, num_last_entries, n);
  }

  assert(idata == num_data); /* index now points to one beyond the end of the array */
  if (data_read != num_data) {
    fprintf(stderr, "Expected %d data points but read %d: problem with file.\n", num_data, data_read);
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

/* detect protein-occupied space by dielectric constant a bit greater than PDIE
   (could also check for SDIE...)
 */
#define PDIE_FUDGE_DELTA 0.05

float draw_diel(const t_membrane *M, float *diel, const float x, const float y, const float z) {
  float R, R_temp;
  R = sqrt(SQR(x - M->x0_p) + SQR(y - M->y0_p));	
  R_temp = (M->R_m1*(z - M->z_m0) - M->R_m0*(z - M->z_m1))/(M->z_m1 - M->z_m0);  	

  if (z >= M->z_m0 && z <= M->z_m1 && *diel > M->pdie+PDIE_FUDGE_DELTA) {
    /* inside membrane region but NOT where the protein is */
    /* bilayer or headgroup dielectric constant outside channel */
    *diel = (z >= M->z_h0 && z <= M->z_h1) ? M->mdie : M->idie; 
    if (R <= R_temp) {
      *diel = M->cdie;   /* channel dielectric painted over membrane */
    }
  }
  return *diel;
}

void print_help() {
  printf("\nusage: draw_membrane2a [OPTIONS] infix\n"
	 "\n"
	 "This program takes the dielectric, kappa, and charge maps from APBS and\n"
	 "accounts for the membrane. The thickness (-d) and the bottom of the\n"
	 "membrane (-z) must be given. We assume that the membrane normal is\n"
	 "along the z-axis. This requires lining the protein along the z-axis\n"
	 "before running this program.  We also output a charge_change_map that\n"
	 "tells me which positions in the charge matrix were edited by the\n"
	 "addition of the membrane.\n"
	 "If a headgroup region (-a l_h -i IDIE) is selected, then the hydrophobic core\n"
	 "of the membrane (eps=MDIE) is modelled with a thickness of l_m - 2*l_h.\n"
	 "\n"
	 "ARGUMENTS:\n"
	 "\n"
	 "  infix  - all maps are constructed as <name><infix>.dx[.gz] where <name> is\n"
	 "           hard-coded (dielx,diely,dielz,kappa,charge); see also OUTPUTS.\n"
	 "\n"
	 "OPTIONS:\n"
	 "  -h          this help\n"
	 "  -z z_m0     bottom of the membrane [-20]\n"
	 "  -d l_m      total membrane thickness (Angstrom) [40]\n"
	 "  -a l_h      headgroup thickness (Angstrom) [0]\n"
	 "  -V V        cytoplasmic potential (kT/e)  [0] UNTESTED\n"
	 "  -I I        molar conc. of one salt-species [0.1]\n"
	 "  -R R_m1     excl. radius at top of  membrane [0]\n" 
	 "  -r R_m0     excl. radius at  bottom of membrane [0]\n"
	 "  -p PDIE     protein dielectric constant (must be SAME as in APBS!) [10]\n"
	 "  -s SDIE     solvent dielectric (only used as default for CDIE) [80]\n"   
	 "  -c CDIE     channel dielectric for (z_m0,R_m0)->(z_m0+l_m,R_m1) [SDIE]\n"   
	 "  -m MDIE     membrane dielectric [2]\n"            
	 "  -i IDIE     headgroup dielectric [MDIR]\n"            
	 "  -Z          read and write gzipped files\n"                    
	 "\n"	                                               
	 "OUTPUTS:\n"
	 "  maps - names are <name><infix>m.dx[.gz]\n"
	 "\n");
}

void print_membrane(const t_membrane *M) {
  if (M->idie != M->mdie && M->z_m0 != M->z_h0) {
    /* we model headgroups */
    printf("Membrane geometry (with headgroups):\n"
	   "\n"
	   "  -------------------- z_m1 = %.1f\n"
	   "     headgroups        idie = %.1f\n"
	   "  ==================== z_h1 = %.1f\n"
	   "     hydrophobic       mdie = %.1f\n"
	   "        core                      \n"
	   "  ==================== z_h0 = %.1f\n"
	   "     headgroups        idie = %.1f\n"
	   "  -------------------- z_m0 = %.1f\n"
	   "\n",
	   M->z_m1, M->idie, M->z_h1, M->mdie, M->z_h0, M->idie, M->z_m0);
  }
  else {
    /* just the slab */
    printf("Membrane geometry (slab only):\n"
	   "\n"
	   "  ==================== z_m1 = %.1f\n"
	   "     hydrophobic       mdie = %.1f\n"
	   "        core                      \n"
	   "  ==================== z_m0 = %.1f\n"
	   "\n",
	   M->z_m1, M->mdie, M->z_m0);
  }
}

void print_exclusionzone(const t_membrane *M) {
  if (M->cdie != M->mdie && (M->R_m1 != 0 || M->R_m0 != 0)) {
    printf("Channel exclusion zone:\n\n");
    if (M->R_m1 > M->R_m0) {
      printf("   ******     R_top = %.1f A   z_m1 = %.1f\n"
	     "    ****      cdie  = %.1f\n"
	     "     **       R_bot = %.1f A   z_m0 = %.1f\n",
	     M->R_m1, M->z_m1, M->cdie, M->R_m0, M->z_m0);
    }
    else if (M->R_m1 < M->R_m0) {
      printf("     **       R_top = %.1f A   z_m1 = %.1f\n"
	     "    ****      cdie  = %.1f\n"
	     "   ******     R_bot = %.1f A   z_m0 = %.1f\n",
	       M->R_m1, M->z_m1, M->cdie, M->R_m0, M->z_m0);
    }
    else {
      printf("     ***      R_top = %.1f A   z_m1 = %.1f\n"
	     "     ***      cdie  = %.1f\n"
	     "     ***      R_bot = %.1f A   z_m0 = %.1f\n",
	     M->R_m1, M->z_m1, M->cdie, M->R_m0, M->z_m0);
    }
    printf("\n");
  }
  else {
    printf("No channel exclusion zone defined.\n");
  }
}

/* access argv after getopt
   optind is a global
 */
int get_argv_int(char **argv, const int index) {
  return atoi(argv[index + optind]);
}

int get_argv_float(char **argv, const int index) {
  return atof(argv[index + optind]);
}

char *get_argv_str(char **argv, const int index) {
  return argv[index + optind];
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
  float tmp_x,tmp_y,tmp_z,dx,dy,dz,l_c_x,l_c_y,l_c_z;
  float x0_p, y0_p, z0_p;
  float x0_x, y0_x, z0_x;
  float x0_y, y0_y, z0_y; 
  float x0_z, y0_z, z0_z; 
  float x0, y0, z0;
  float V=0, I=0.1, sdie, cdie, pdie, mdie, idie;
  float l_m=40, a0=0;
  float z_m0=0, R_m0=0, R_m1=0;
  float R, R_temp;
  char infix[MAXLEN];
  char *file_name_x, *file_name_y, *file_name_z;
  char *file_name_k, *file_name_c;
  char *f1, *f2, *f3, *f4, *f5, *f6;
  char ext[]="m.dx";
  bool compression = FALSE;
  int c;
  t_membrane Membrane;
 

  printf("----------------------------------------------------------------\n");
  printf("* draw_membrane2a.c                                   12/01/10 *\n");  /* magic version line */
  printf("----------------------------------------------------------------\n");
  printf("draw_membrane2  -- (c) 2008 Michael Grabe [09/02/08]\n");
  printf("draw_membrane4  -- (c) 2010 Michael Grabe [07/29/10]\n");
  printf("draw_membrane2a -- (c) 2010 Oliver Beckstein (options&gzipped files) [12/01/10]\n");
  printf("Published under the Open Source MIT License (see http://sourceforge.net/projects/apbsmem/).\n");
  printf("Based on http://www.poissonboltzmann.org/apbs/examples/potentials-of-mean-force/the-polar-solvation-potential-of-mean-force-for-a-helix-in-a-dielectric-slab-membrane/draw_membrane2.c\n");
  printf("----------------------------------------------------------------\n");
  
  /* explicit defaults for options */
  mdie = 2.0;    /* watch out for this used to be 10.0 */ 
  idie = -1;     /* set to mdie iff < 0 */
  sdie = 80.0;
  cdie = -1;     /* set to sdie iff < 0 */
  pdie = 10.0;   /* MUST match the value set in APBS !! */
  
  z_m0 = -20;    /* lower z of membrane */  
  l_m = 40;      /* thickness of membrane */

  R_m1 = 0;      /* upper exclusion radius */
  R_m0 = 0;      /* lower exclusion radius */

  opterr = 0;
  while ((c = getopt(argc, argv, "hZz:d:s:c:m:p:i:a:V:I:r:R:")) != -1) {
    switch(c) {
    case 'h':
      print_help();
      return EXIT_SUCCESS;
    case 'z':
      z_m0 = atof(optarg);
      break;
    case 'd':
      l_m = atof(optarg);
      break;
    case 'a':
      a0 = atof(optarg);
      break;
    case 's':
      sdie = atof(optarg);
      break;
    case 'c':
      cdie = atof(optarg);
      break;
    case 'm':
      mdie = atof(optarg);
      break;
    case 'i':
      idie = atof(optarg);
      break;
    case 'p':
      pdie = atof(optarg);
      break;
    case 'V':
      V = atof(optarg);  /* membrane potential in kT/e--- TODO: use mV as input */
      break;
    case 'I':
      I = atof(optarg);  /* ionic strength in mol/l */
      break;
    case 'r':
      R_m0 = atof(optarg);
      break;
    case 'R':
      R_m1 = atof(optarg);
      break;     
    case 'Z':
      compression = TRUE;
      break;
    case '?':
      if (optopt == 'z' || optopt == 'd' || optopt == 's' || optopt == 'c' ||
	  optopt == 'i' || optopt == 'a' ||
	  optopt == 'm' || optopt == 'p' || optopt == 'V' || optopt == 'I' || 
	  optopt == 'r' || optopt == 'R')
	fprintf(stderr, "Option -%c requires an argument.\n", optopt);
      else 
	fprintf(stderr, "Unknown option `-%c'.\n", optopt);  // should check isprint(optopt)...
      return EXIT_FAILURE;
    default:
      abort();
    }
  }

  if (argc-optind != 1) {
    fprintf(stderr, "%d argument%s were provided on the commandline "
	    "but exactly one (the infix) is required. See `-h' for help.\n", 
	    argc-optind, argc-optind==1 ? "":"s");
    return EXIT_FAILURE;
  }

  strcpy(infix, get_argv_str(argv, 0));
  printf("Using hard-coded names with your infix to find files: infix=%s\n", infix);

  if (compression)
    printf("Reading gzip-compressed dx files.\n");
  else
    printf("Reading un-compressed dx files (default).\n");

  if (cdie < 0)
    cdie = sdie;
  if (cdie != sdie && (R_m0 > 0 || R_m1 > 0))
    printf("Using different dielectric in channel (%.1f) than in bulk (%.1f)\n", cdie, sdie);
  
  if (idie < 0)
    idie = mdie;
  if (idie != mdie || a0 > 0)
    printf("Modelling headgroup region of thickness %.1f with epsilon=%.1f\n", a0, idie);

  printf("Running with these arguments:\n");
  printf(">>> draw_membrane2a %s -s %.1f -c %.1f -m %.1f -p %.1f -i %.1f -a %.1f "
	 "-V %.3f -I %.3f -z %.3f -d %.3f -r %.1f -R %.1f  %s\n",
	 compression ? "-Z" : "", sdie, cdie, mdie, pdie, idie, a0, V, I,
	 z_m0, l_m, R_m0, R_m1, infix);

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

  /*****************************************************/
  /* read in the x-shifted dielectric data             */
  /*****************************************************/

  in = gzopen(file_name_x,"r");
  if (in == NULL) {
    fprintf(stderr, "ERROR: %s not found in current directory.\n",  file_name_x);
    fprintf(stderr, "       See -h for usage.\n");
    return EXIT_FAILURE;
  }
  printf("Reading %s...\n", file_name_x);

  /* First read the header */
  read_header(in, &dim_x, &dim_y, &dim_z, &x0_x, &y0_x, &z0_x, &dx, &dy, &dz, &dim3);

  /* assign the memory to the arrays */
  x_x = calloc(dim_x+1,sizeof(float));
  y_x = calloc(dim_y+1,sizeof(float));
  z_x = calloc(dim_z+1,sizeof(float));
  x_y = calloc(dim_x+1,sizeof(float));
  y_y = calloc(dim_y+1,sizeof(float));
  z_y = calloc(dim_z+1,sizeof(float));
  x_z = calloc(dim_x+1,sizeof(float));
  y_z = calloc(dim_y+1,sizeof(float));
  z_z = calloc(dim_z+1,sizeof(float));
  d_x = calloc(dim_x*dim_y*dim_z+1,sizeof(float));
  d_y = calloc(dim_x*dim_y*dim_z+1,sizeof(float));
  d_z = calloc(dim_x*dim_y*dim_z+1,sizeof(float));
  /* Now the Kappa and charge Arrays */
  x = calloc(dim_x+1,sizeof(float));
  y = calloc(dim_y+1,sizeof(float));
  z = calloc(dim_z+1,sizeof(float));
  kk = calloc(dim_x*dim_y*dim_z+1,sizeof(float));
  cc = calloc(dim_x*dim_y*dim_z+1,sizeof(float));
  map = calloc(dim_x*dim_y*dim_z+1,sizeof(int));

  /* initialize x,y,z, and diel vectors */
  l_c_x = dim_x*dx;
  l_c_y = dim_y*dx;
  l_c_z = dim_z*dx;

  tmp_x = x0_x;
  tmp_y = y0_x;
  tmp_z = z0_x;

  for (i=1; i <= dim_x; ++i) {
    x_x[i]=tmp_x;
    tmp_x+=dx;
  }
  for (i=1; i <= dim_y; ++i) {
    y_x[i]=tmp_y;
    tmp_y+=dy;
  }
  for (i=1; i <= dim_z; ++i) {
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
  x0_p = x0_x+l_c_x/2-dx/2;  /* this is the shift term that moves half-step off grid */ 
  y0_p = y0_x+l_c_y/2;
  z0_p = z0_x+l_c_z/2;

  /* TODO: use t_membrane throughout, currently only used for draw() */
  Membrane.z_m0 = z_m0;         /* bottom of membrane */
  Membrane.z_m1 = z_m0 + l_m;   /* top of the membrane */
  Membrane.z_h0 = z_m0 + a0;    /* top of lower head group region */
  Membrane.z_h1 = Membrane.z_m1 - a0;    /* bottom of upper head group region */
  Membrane.pdie = pdie;
  Membrane.mdie = mdie;
  Membrane.cdie = cdie;
  Membrane.idie = idie;
  Membrane.R_m0 = R_m0;
  Membrane.R_m1 = R_m1;
  Membrane.x0_p = x0_p;
  Membrane.y0_p = y0_p;
  Membrane.z0_p = z0_p;

  print_membrane(&Membrane);
  print_exclusionzone(&Membrane);

  /*****************************************************/
  /* read in the y-shifted dielectric data             */
  /*****************************************************/

  in = gzopen(file_name_y,"r");
  if (in == NULL) {
    printf("File name %s not found.\n", file_name_y);
    return EXIT_FAILURE;
  }
  printf("Reading %s...\n", file_name_y);

  /* First read the header */
  read_header(in, &dim_x, &dim_y, &dim_z, &x0_y, &y0_y, &z0_y, &dx, &dy, &dz, &dim3);

  /* initialize x,y,z, and diel vectors */
  tmp_x = x0_y;
  tmp_y = y0_y;
  tmp_z = z0_y;

  for (i=1; i <= dim_x; ++i) {
    x_y[i]=tmp_x;
    tmp_x+=dx;
  }
  for (i=1; i <= dim_y; ++i) {
    y_y[i]=tmp_y;
    tmp_y+=dy;
  }
  for (i=1; i <= dim_z; ++i) {
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
    return EXIT_FAILURE;
  }
  printf("Reading %s...\n", file_name_z);

  /* First read the header */
  read_header(in, &dim_x, &dim_y, &dim_z, &x0_z, &y0_z, &z0_z, &dx, &dy, &dz, &dim3);

  /* initialize x,y,z, and diel vectors */
  tmp_x = x0_z;
  tmp_y = y0_z;
  tmp_z = z0_z;

  for (i=1; i <= dim_x; ++i) {
    x_z[i]=tmp_x;
    tmp_x+=dx;
  }
  for (i=1; i <= dim_y; ++i) {
    y_z[i]=tmp_y;
    tmp_y+=dy;
  }
  for (i=1; i <= dim_z; ++i) {
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
    return EXIT_FAILURE;
  }
  printf("Reading %s...\n", file_name_k);

  /* First read the header */
  read_header(in, &dim_x, &dim_y, &dim_z, &x0, &y0, &z0, &dx, &dy, &dz, &dim3);

  /* initialize x,y,z, and kappa vectors */
  tmp_x = x0;
  tmp_y = y0;
  tmp_z = z0;

  for (i=1; i <= dim_x; ++i) {
    x[i]=tmp_x;
    tmp_x+=dx;
  }
  for (i=1; i <= dim_y; ++i) {
    y[i]=tmp_y;
    tmp_y+=dy;
  }
  for (i=1; i <= dim_z; ++i) {
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
    return EXIT_FAILURE;
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
  printf("Adding the membrane to the maps...\n");

  /******************************************************/
  /* set up the vector                                  */
  /******************************************************/
  /* NOTE: Array indexing is done by FORTRAN convention (1-based, z fastest varying) */
  cnt=1;

  for (k=1; k <= dim_x; ++k) {  /* loop over x */
    for (j=1; j <= dim_y; ++j) {  /* loop over y */
      for (i=1; i <= dim_z; ++i) {  /* loop over z */
	draw_diel(&Membrane, &(d_x[cnt]), x_x[k], y_x[j], z_x[i]);
	draw_diel(&Membrane, &(d_y[cnt]), x_y[k], y_y[j], z_y[i]);
	draw_diel(&Membrane, &(d_z[cnt]), x_z[k], y_z[j], z_z[i]);

	if (z[i] <= z_m0 && kk[cnt] != 0.0) {
	  /* charge for mem V */
	  /* see my notes for this expression */
	  cc[cnt] = 0.0012045*I*V;
	  /* update the change map */
	  map[cnt] = 1;
	}
	else {
	  /* position was not changed */
	  map[cnt] = 0;
	}

	R = sqrt(SQR(x[k]-Membrane.x0_p) + SQR(y[j]-Membrane.y0_p));                      
	R_temp = (R_m1*(z[i]-Membrane.z_m0) - R_m0*(z[i]-Membrane.z_m1))/(Membrane.z_m1 - Membrane.z_m0);

	if (z[i] <= Membrane.z_m1 && z[i] >= Membrane.z_m0 && R > R_temp) {
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

  /********************X-DATA******************************/
  /* add the "m" extension to the file */
  f1 = newname("dielx", infix, ext, compression);
  out = xopen(compression, f1,"w");
  write_header(compression, out, "X-SHIFTED DIELECTRIC MAP", z_m0, l_m, 
	       dim_x, dim_y, dim_z, x0_x, y0_x, z0_x, dx ,dy, dz);
  write_data(compression, out, dim3, d_x);
  write_attr_positions(compression, out);
  xclose(compression, out);
  printf("Wrote %s.\n", f1);

  /********************Y-DATA******************************/
  /* give the file an "m" extension */
  f2 = newname("diely", infix, ext, compression);
  out = xopen(compression, f2,"w");
  write_header(compression, out, "Y-SHIFTED DIELECTRIC MAP", z_m0, l_m, 
	       dim_x, dim_y, dim_z, x0_y, y0_y, z0_y, dx ,dy, dz);
  write_data(compression, out, dim3, d_y);
  write_attr_positions(compression, out);
  xclose(compression, out);
  printf("Wrote %s.\n", f2);

  /**********************Z-DATA*****************************/
  /* give the file an "m" extension */
  f3 = newname("dielz", infix, ext, compression);
  out = xopen(compression, f3,"w");
  write_header(compression, out, "Z-SHIFTED DIELECTRIC MAP", z_m0, l_m, 
	       dim_x, dim_y, dim_z, x0_z, y0_z, z0_z, dx ,dy, dz);
  write_data(compression, out, dim3, d_z);
  write_attr_positions(compression, out);
  xclose(compression, out);
  printf("Wrote %s.\n", f3);

  /*********************KAPPA******************************/
  /* give the file an "m" extension */
  f4 = newname("kappa", infix, ext, compression);
  out = xopen(compression, f4,"w");
  write_header(compression, out, "KAPPA MAP", z_m0, l_m, dim_x, dim_y, dim_z, x0, y0, z0, dx ,dy, dz);
  write_data(compression, out, dim3, kk);
  write_attr_positions(compression, out);
  xclose(compression, out);
  printf("Wrote %s.\n", f4);

  /********************CHARGE*******************************/
  /* give the file an "m" extension */
  f5 = newname("charge", infix, ext, compression);
  out = xopen(compression, f5,"w");
  write_header(compression, out, "CHARGE MAP", z_m0, l_m, dim_x, dim_y, dim_z, x0, y0, z0, dx ,dy, dz);
  write_data(compression, out, dim3, cc);
  write_attr_positions(compression, out);
  xclose(compression, out);
  printf("Wrote %s.\n", f5);

  /********************CHARGE CHANGE MAP*************************/
  f6 = newname("change_map", infix, ext, compression);
  out = xopen(compression, f6,"w");
  write_header(compression, out, "CHARGE CHANGE MAP", z_m0, l_m, dim_x, dim_y, dim_z, x0, y0, z0, dx ,dy, dz);
  write_data(compression, out, dim3, map);  // int map is converted to float
  write_attr_positions(compression, out);
  xclose(compression, out);
  printf("Wrote %s.\n", f6);

  /***********************************************************/
  /* clean up */
  free(x_x); free(y_x); free(z_x);
  free(x_y); free(y_y); free(z_y); 
  free(x_z); free(y_z); free(z_z); 
  free(d_x); free(d_y); free(d_z);
  free(x);   free(y);   free(z);
  free(kk);  free(cc);  free(map);

  free(f1); free(f2); free(f3); free(f4); free(f5); free(f6); /* calloc'ed by newname() :-p */

  printf("Your files have been written.\n");
  return EXIT_SUCCESS;
}
