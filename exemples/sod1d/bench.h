/* simple benchmark header only

-- USAGE : 
    `BENCH_START` : start the timer
    `BENCH_LOG`   : log the time, increase iteration counter
    `BENCH_MEAN(*m)` : get mean
    `BENCH_MEAN_STD(*m,*s)` : get mean and stdev
    `BENCH_LAST`  : get last logged time
    `BENCH_PRINT` : print mean and stdev
    `BENCH_PRINT_PERIOD(t)` : print every 't' iterations
    `BENCH_PERIOD(t)` : evaluate to true 't' iterations
    `BENCH_BENCH(nb_repeat, function)` : bench `function` with '*nb_repeat*' repetition
    `BENCH_BENCH_PRINT(nb_repeat, print_period, function)` : bench `function` and print every '*print_period*' iterations
    `BENCH_RESET`     : reset iteration counter
    `BENCH_INIT` : init necessary variables

    warning : mean and stdev are incorrect until the array is filled once (a star* is put at end of the print to warn of this case).

-- CONFIG :
    `BENCH_LOG_SIZE`  : set the log size (default : 128)
    `BENCH_PRECISION` : set the time precision : s, ms, µs, ns (default : ms)
    `BENCH_fd` : set the output FILE* of fprintf (default : stdout)
    `BENCH_PRINT_FORMAT` : printf format for time (default : "%9.4lf")
    `BENCH_NO_WARN_LOGSIZE` : disable the star* at end of the print to warn that the array isn't filled.
    `BENCH_CSV_PRINT` : set the print output in csv format (mean, stdev, iterations)
    `BENCH_NO_AUTO_INIT` : disable init globals variables

-- DEPENDENCIES :
    `math.h` : sqrt(), define `BENCH_NO_SQRT` to disable and get variance 
               instead of stdev, or define your own sqrt in `__BENCH_SQRT__`
    `time.h` : clock_gettime()
*/

#include <time.h>

#ifndef BENCH_LOG_SIZE
#define BENCH_LOG_SIZE 128
#endif
#ifndef BENCH_PRECISION
#define BENCH_PRECISION ms
#endif
#ifndef BENCH_fd
#define BENCH_fd stdout
#endif
#ifndef BENCH_PRINT_FORMAT
#define BENCH_PRINT_FORMAT "%9.4lf"
#endif
#ifndef BENCH_NO_WARN_LOGSIZE
#define __BENCH_WARN_LOGSIZE__                                 \
            (__bench_itetation__ < BENCH_LOG_SIZE) ? "*" : ""
#else
#define __BENCH_WARN_LOGSIZE__  ""
#endif
#ifndef BENCH_CSV_PRINT
#define BENCH_PRINT __BENCH_PRINT__
#else
#define BENCH_PRINT __BENCH_CPRINT__
#endif


#if !defined(BENCH_NO_SQRT) && !defined(__BENCH_SQRT__)
#include <math.h>
#define __BENCH_SQRT__ sqrt
#elif !defined(__BENCH_SQRT__)
#define __BENCH_SQRT__
#endif

#define __BENCH_s__     1.
#define __BENCH_ms__    1e3
#define __BENCH_us__    1e6
#define __BENCH_ns__    1e9
#define __BENCH_µs__    __BENCH_us__

#define BENCH_INIT                                                 \
            __attribute_maybe_unused__ unsigned long __bench_itetation__ = 0;          \
            __attribute_maybe_unused__ double __bench_times__[BENCH_LOG_SIZE];         \
            __attribute_maybe_unused__ struct timespec __bench_start__, __bench_end__;

#define BENCH_START \
            clock_gettime(CLOCK_MONOTONIC_RAW, &__bench_start__);

#define BENCH_LOG                                                                               \
        {                                                                                       \
            clock_gettime(CLOCK_MONOTONIC_RAW, &__bench_end__);                                 \
            double __sec__ = (double)(__bench_end__.tv_sec  - __bench_start__.tv_sec)           \
                           + (double)(__bench_end__.tv_nsec - __bench_start__.tv_nsec) * 1e-9;  \
            __bench_times__[__bench_itetation__ % BENCH_LOG_SIZE] = __sec__ * __BENCH_PREC__;   \
            __bench_itetation__++;                                                              \
        }

#define BENCH_BENCH(nb_repeat, function)                        \
        {                                                       \
            BENCH_INIT                                          \
            for(unsigned long i = 0; i < nb_repeat; i++){       \
                BENCH_START                                     \
                function;                                       \
                BENCH_LOG                                       \
            }                                                   \
        }

#define BENCH_BENCH_PRINT(nb_repeat, print_period, function)    \
        {                                                       \
            BENCH_INIT                                          \
            for(unsigned long i = 0; i < nb_repeat; i++){       \
                BENCH_START                                     \
                function;                                       \
                BENCH_LOG                                       \
                BENCH_PRINT_PERIOD(print_period)                \
            }                                                   \
        }

#define BENCH_MEAN(m)           \
        {                       \
            __BENCH_COMPUTE__;  \
            *m = __mean__;      \
        }

#define BENCH_MEAN_STD(m,s)     \
        {                       \
            __BENCH_COMPUTE__;  \
            *m = __mean__;      \
            *s = __std__;       \
        }

#define BENCH_LAST  \
            __bench_times__[(__bench_itetation__-1) % BENCH_LOG_SIZE]

#define __BENCH_PRINT__                                         \
        {                                                       \
            __BENCH_COMPUTE__;                                  \
            fprintf(BENCH_fd,                                   \
                    "%s" BENCH_PRINT_FORMAT                     \
                    " " __BENCH_STR__(BENCH_PRECISION)          \
                    " (± "BENCH_PRINT_FORMAT")%15s %lu %4s\n",  \
                    "bench time : ", __mean__, __std__,         \
                    "iteration :", __bench_itetation__,         \
                    __BENCH_WARN_LOGSIZE__);                    \
        }

#define __BENCH_CPRINT__                       \
        {                                      \
            __BENCH_COMPUTE__;                 \
            fprintf(BENCH_fd,                  \
                    BENCH_PRINT_FORMAT "\t"    \
                    BENCH_PRINT_FORMAT "\t"    \
                    "%lu\t%s\n",               \
                    __mean__, __std__,         \
                    __bench_itetation__,       \
                    __BENCH_WARN_LOGSIZE__);   \
        }

#define BENCH_PRINT_PERIOD(t) \
            if(BENCH_PERIOD(t)) BENCH_PRINT;

#define BENCH_PERIOD(t) \
            (__bench_itetation__ % t == 0)

#define BENCH_RESET \
            __bench_itetation__ = 0;

#ifndef BENCH_NO_AUTO_INIT
    BENCH_INIT
#endif

#define __BENCH_COMPUTE__                                             \
            double __mean__ = 0.;                                     \
            double __std__ = 0.;                                      \
            for(unsigned long i = 0; i < BENCH_LOG_SIZE; i++){        \
                __mean__ += __bench_times__[i];                       \
                __std__  += __bench_times__[i] * __bench_times__[i];  \
            }                                                         \
            __mean__ /= BENCH_LOG_SIZE;                               \
            __std__ /= BENCH_LOG_SIZE;                                \
            __std__ = __BENCH_SQRT__(__std__ - __mean__*__mean__);

#define __BENCH_CAT_HELPER__(x, y, z) x##y##z
#define __BENCH_CAT__(x, y, z) __BENCH_CAT_HELPER__(x, y, z)
#define __BENCH_STR_HELPER__(x) #x
#define __BENCH_STR__(x) __BENCH_STR_HELPER__(x)
#define __BENCH_PREC__ __BENCH_CAT__(__BENCH_, BENCH_PRECISION, __)