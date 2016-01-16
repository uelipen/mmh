/* file fpe_handler.c: trap floating point signals on IRIX 6
                       from consult@ncsa, March 1995 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/fpu.h>
#include <signal.h>
void fpe_abort( void );
#define INEXACT_EXCEPTION       0
#define UNDERFLOW_EXCEPTION     1
#define OVERFLOW_EXCEPTION      2
#define DIV0_EXCEPTION          3
#define INVALID_EXCEPTION       4
#define UNKNOWN_EXCEPTION       5
char *fpe_cause[] = {
        "Inexact Exception",
        "Underflow Exception",
        "Overflow Exception",
        "Division by Zero",
        "Invalid Operation",
        "Unknown Exception"
};
/*  Initialize the "enable bits"  */
#define INEXACT_EXCEPTION_ABORT         1
#define UNDERFLOW_EXCEPTION_ABORT       0
#define OVERFLOW_EXCEPTION_ABORT        1
#define DIV0_EXCEPTION_ABORT            1
#define INVALID_EXCEPTION_ABORT         1
void set_fpe_( ){
    set_fpe( );
}
check_fpe( ) {
    union fpc_csr       f;
    int fpe_exception;
    f.fc_word   = get_fpc_csr();
/* Clear the exception flags */
    if( f.fc_struct.se_divide0 == 1) {
        fpe_exception = DIV0_EXCEPTION ;
    } else if( f.fc_struct.se_invalid == 1) {
        fpe_exception = INVALID_EXCEPTION ;
    } else if( f.fc_struct.se_underflow == 1) {
        fpe_exception = UNDERFLOW_EXCEPTION ;
    } else if( f.fc_struct.se_overflow == 1) {
        fpe_exception = OVERFLOW_EXCEPTION;
    } else if( f.fc_struct.se_inexact == 1) {
        fpe_exception = INEXACT_EXCEPTION;
    } else{
        fpe_exception = UNKNOWN_EXCEPTION;
    }
    return( fpe_exception );
}
set_fpe( ) {
    union fpc_csr       f;
    f.fc_word   = get_fpc_csr();
/* Set up the Enable bits */
/*    f.fc_struct.en_inexact   = INEXACT_EXCEPTION_ABORT ;
    f.fc_struct.en_underflow = UNDERFLOW_EXCEPTION_ABORT ;
 */
    f.fc_struct.en_overflow  = OVERFLOW_EXCEPTION_ABORT ;
    f.fc_struct.en_divide0   = DIV0_EXCEPTION_ABORT ;
    f.fc_struct.en_invalid   = INVALID_EXCEPTION_ABORT ;
/* Clear the exception flags */
    f.fc_struct.se_inexact   = 0 ;
    f.fc_struct.se_underflow = 0 ;
    f.fc_struct.se_overflow  = 0 ;
    f.fc_struct.se_divide0   = 0 ;
    f.fc_struct.se_invalid   = 0 ;
    set_fpc_csr(f.fc_word);
    signal( SIGFPE, fpe_abort );
    return( f.fc_word );
}
/*  Kill the Process with a SIGFPE signal */
void fpe_abort( void ){
    fprintf(stderr,
        "Floating Point Exception of type: %s -ABORTING\n",
        fpe_cause[check_fpe()]);
    kill( getpid(), SIGFPE );
}


#define MAX(x,y) ((x>y)?x:y)
void csysload_(double *sysload)
{
    char *cmd = "/usr/bsd/uptime | awk '{print $(NF-2),$(NF-1),$(NF)}'";
    FILE *ptr;
    int errcode;
    double ld1,ld5,ld15;
 
    if ((ptr = popen(cmd, "r")) != NULL) {
        fscanf(ptr,"%lf, %lf, %lf",&ld1,&ld5,&ld15);
        *sysload=MAX(ld1,MAX(ld5,ld15));
        errcode=pclose(ptr);
        if (errcode != 0)
           fprintf(stderr,"csysload: %s returned %d\n",cmd,errcode);
    } else {
        fprintf(stderr,"csysload: error reading system command\n");
        *sysload=-1.;
    }
}

void cncpu_(int *ncpu)
{
    char *cmd = "hinv -c processor | grep 'MHZ' | wc -l";
    /*  "hinv -c processor | head -1 | awk '{print $1}'"; */
    FILE *ptr;
    int errcode, ntcpu;
    double ld1,ld5,ld15;
 
    if ((ptr = popen(cmd, "r")) != NULL) {
        fscanf(ptr,"%d",&ntcpu);
        errcode=pclose(ptr);
        if (errcode != 0)
           fprintf(stderr,"cncpu: %s returned %d\n",cmd,errcode);
        else if (ntcpu > 1 && ntcpu < 100) *ncpu=ntcpu;
    } else {
        fprintf(stderr,"cncpu: error reading system command\n");
    }
}
