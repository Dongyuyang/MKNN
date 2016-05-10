#ifndef _CATCH_H_
#define _CATCH_H_

/*
 * Function: Definition for catching run-time in processing.
 * Author  : Jianquan LIU, DSE Lab., Univ. of Tsukuba, Japan.
 * Data    : Tue May 29 20:36:21 JST 2007
 * Email   : ljq@dblab.is.tsukuba.ac.jp
 */

#include <iostream>
using namespace std;

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <sys/time.h>

#define DEFAULT_DESC "time cost"
#define UNIT_MIN     "minute(s)"
#define UNIT_SEC     "second(s)"
#define UNIT_MSEC    "millisecond(s)"
#define UNIT_USEC    "microsecond(s)"

#define USEC_PER_MIN  (60000000.0)
#define USEC_PER_SEC   (1000000.0)
#define USEC_PER_MSEC     (1000.0)

#ifndef STRLEN
#define STRLEN        (256)
#endif



class CATCH{

    public:
        CATCH(){
           timerclear(&m_old);
           timerclear(&m_cur);
        }
        ~CATCH(){
        }

        //definition for the time-type of displaying
        enum {
            IN_MIN=0,  //display time cost in unit of minute
            IN_SEC=1,  //display time cost in unit of second
            IN_MSEC=2, //display time cost in unit of millisecond
            IN_USEC=3  //display time cost in unit of microsecond
        };

        //menber functions
        void catch_time();
        double get_cost(const int t=IN_SEC);
        void print_time(const char *desc=DEFAULT_DESC,  //description for display
                        const int t=IN_SEC);             //time-type for display

    private:
        struct timeval m_old;
        struct timeval m_cur;

};

void CATCH::catch_time()
{
    if (!timerisset(&m_old))
    {
        gettimeofday(&m_old, NULL);
        m_cur = m_old;
    }
    else
    {
        m_old = m_cur;
        gettimeofday(&m_cur, NULL);
    }

    return ;
}

double CATCH::get_cost(const int t)
{
    double n_sec   = (double)(m_cur.tv_sec - m_old.tv_sec);
    double n_usec  = (double)(m_cur.tv_usec - m_old.tv_usec);
    double n_total = n_sec * 1000000.0 + n_usec;

    switch(t)
    {
        case IN_MIN:
            n_total /= USEC_PER_MIN;
            break;

        case IN_MSEC:
            n_total /= USEC_PER_MSEC;
            break;

        case IN_USEC:
            break;

        case IN_SEC:
        default    :
            n_total /= USEC_PER_SEC;
            break;
    }

    return n_total;
}

void CATCH::print_time(const char *desc, const int t)
{
    double n_sec   = (double)(m_cur.tv_sec - m_old.tv_sec);
    double n_usec  = (double)(m_cur.tv_usec - m_old.tv_usec);
    double n_total = n_sec * 1000000.0 + n_usec;
    char sformat[STRLEN] = {0};

    switch(t)
    {
        case IN_MIN:
            snprintf(sformat, STRLEN, "%s : %%.2f %s\n\n", desc, UNIT_MIN);
            n_total /= USEC_PER_MIN;
            break;

        case IN_MSEC:
            snprintf(sformat, STRLEN, "%s : %%.2f %s\n\n", desc, UNIT_MSEC);
            n_total /= USEC_PER_MSEC;
            break;

        case IN_USEC:
            snprintf(sformat, STRLEN, "%s : %%.2f %s\n\n", desc, UNIT_USEC);
            break;

        case IN_SEC:
        default    :
            snprintf(sformat, STRLEN, "%s : %%.2f %s\n\n", desc, UNIT_SEC);
            n_total /= USEC_PER_SEC;
            break;
    }

    printf(sformat, n_total);

    return ;
}

#endif
