#include <stdio.h>
#include <stdlib.h>
#include "aceplug.h"
#include "metadyn.h"

int aceplug_init( 
	struct aceplug_sim_t *s,
	int argc,
	char **argkey,
	char **argval)
{   
  int i;
  int nrepl,repl;
  real dt;
  real rte0,rteio;
  real *mass;
  real *charge;
  char *metainp=NULL;
  real box[3];

  if(s->version < 200 || s->version > 300) {
    printf("PLUMED: wrong ACEMD plugin interface version (%d), expected 200 to 300.\nPlease download a recent ACEMD version from www.acellera.com or compile the 'legacy' version of the plugin.\n", s->version );
    exit(1);
  }

  for(i=0; i < argc; i++ ) {
    if(!strcasecmp(argkey[i],"input")) {
      metainp=argval[i];
    } else {
      printf("PLUMED: ignoring key %s value %s\n", argkey[i], argval[i]);
    }
  }

  if(!metainp) {
    printf("PLUMED: mandatory plugin argument 'input' missing\n");
    exit(1);
  }

  dt = s->timestep_fs;		/* Timestep */

  box[0]=s->box.x;		/* Box */
  box[1]=s->box.y;
  box[2]=s->box.z;

  nrepl = 1;
  repl = 1;

  //     Target_temperature controllare che non sia in kcal/mol (*0.59227/298)
  rte0  = 298.0;
  rteio = rte0; 

  /* We maybe could recycle the pointers as well */
  charge    = (real *)calloc(s->natoms,sizeof(real));
  mass      = (real *)calloc(s->natoms,sizeof(real));
  for(i=0;i<s->natoms;i++) {
    mass[i]   =  s->mass[i];
    charge[i] =  s->charge[i];
  }

  printf("PLUMED: initializing with control file '%s', initial box (%f,%f,%f), timestep %f fs\n",metainp,box[0],box[1],box[2],dt);

  init_metadyn(s->natoms, charge, mass, 
	       dt, repl, nrepl, rte0, rteio, metainp, box );
  return 0;
}


int aceplug_update( struct aceplug_sim_t *s) {
  meta_force_calculation(s);
  return 0;
}

int aceplug_terminate( void *privdata ) {
  printf("PLUMED: exiting\n");
  return 0;
}
