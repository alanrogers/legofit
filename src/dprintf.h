#ifndef LDPSIZ_DPRINTF_H
#  define LDPSIZ_DPRINTF_H
#  include <pthread.h>
#  include <stdio.h>
#  include <string.h>
/// If DPRINTF_ON is defined, DPRINTF((arg)) will acquire a mutex lock,
/// hand arg to printf, flush stdout, and release the lock. On error,
/// it prints an message and exits.
///
/// If DPRINTF_ON is not define, DPRINTF is a noop.
///
/// The mutex "outputlock" is defined in dprintf.c, which must be linked
/// into executables that use DPRINTF.
#  ifdef DPRINTF_ON
#    define DPRINTF(arg) do{                                    \
        int  dpr_status;                                        \
        char dpr_buff[50];                                      \
        dpr_status=pthread_mutex_lock(&outputLock);             \
        if(dpr_status) {                                        \
            strerror_r(dpr_status, dpr_buff, sizeof(dpr_buff)); \
            fprintf(stderr,"%s:%s:%d: lock %d (%s)\n",          \
                    __FILE__,__func__,__LINE__,                 \
                    dpr_status, dpr_buff);                      \
            exit(1);                                            \
        }                                                       \
        printf arg ;                                            \
        fflush(stdout);                                         \
        dpr_status = pthread_mutex_unlock(&outputLock);         \
        if(dpr_status) {                                        \
            strerror_r(dpr_status, dpr_buff, sizeof(dpr_buff)); \
            fprintf(stderr,"%s:%s:%d: unlock %d (%s)\n",        \
                    __FILE__,__func__,__LINE__,                 \
                    dpr_status, dpr_buff);                      \
            exit(1);                                            \
        }                                                       \
    }while(0)
#  else
#    define DPRINTF(arg)
#  endif
#endif
