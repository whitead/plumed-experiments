#ifndef ___PLUGIN
#define ___PLUGIN 1

#if !defined(__VECTOR_TYPES_H__) 

typedef __attribute__((aligned(16))) struct {
  double x;
  double y;
  double z;
  double w;
} double4;

typedef __attribute__((aligned(16))) struct {
  float x;
  float y;
  float z;
  float w;
} float4;
#endif


  /* Time unit in femtoseconds */
#define   TIMEFACTOR   48.88821 
  /* Coulomb's constant for electrostatics, units kcal*A/mol/e^2 */
#define MD_COULOMB  332.0636
  /* Boltzman's constant for temperature, units kcal/A/K */
#define MD_BOLTZMAN  0.001987191

#define MD_KCAL_TO_KJ  4.1868


#define CONCAT(x,y,z) x##y##z
#define CONCAT2(x,y,z) CONCAT(x,y,z)
#define _FUNCTION(prefix,suffix) CONCAT2(prefix,_,suffix)


struct aceplug_sim_t {
	int version;
	double4 *pos;
	float4 *frc;
	float  *mass;
	float  *charge;
	int    natoms;
	unsigned long step;
	void *privdata;
	float4 box;
	float timestep_fs;
} ;


#endif

#ifdef __cplusplus
extern "C" {
#endif

int _FUNCTION(aceplug,init)( struct aceplug_sim_t*, int argc, char **argkey, char **argval );

int _FUNCTION(aceplug,update)( struct aceplug_sim_t* );

int _FUNCTION(aceplug,terminate)( void * privdata );

#ifdef __cplusplus
}
#endif



