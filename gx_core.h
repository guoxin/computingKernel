/*
*  This headfile gives some general support on pthread programmes.
*  
*  CORENUM macro defines in general how many threads are used.
*  
*  the sfgd functions checks whether pthreads functions run normally
*  and give warnings (instead of terminating the whole programme).
*  
*  history: first version: simply the CORENUM macro definition.
*  	 second version: Oct 23, 2011, sfgd functions added.
*        third version: Mar 22, 2012, R.h error exception added.
*/
#define CORENUM 24
#include <R.h>
#include <stdio.h>
#include <stdlib.h>
const static char * gx_default_information = "unspecified";
void gx_sfgd_info(const int state, const char* info){
	/* "sfgd" = safeguard, a trademark of soap in China. */
	if(state != 0){
		printf("\n\nA fatal ERROR about pthread is detected.\n");
		printf("Related state code: %d\n", state);
		printf("Related information: %s\n", info);
		printf("The following output might not be reliable!\n\n");
                error("fatal system error from C code by guo xin.");
		return;
	}else return;
}
void gx_sfgd(const int state){
	if(state != 0){
		gx_sfgd_info(state, gx_default_information);
		return;
	}else return;
}

