/***
 *  File: getWinTime.hpp
 *  Developed: 2008
 *
 *  Author: Xiao Yang <isuxyang@gmail.com>
 *  Copyright (c) 2008 Xiao Yang
 *  Distributed under the Boost Software License.
 *  See accompanying file LICENSE.
 */

#if !defined(GETWINTIME_HPP)
#define GETWINTIME_HPP

#include <time.h>
#include <windows.h>

#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif

struct timezone {
    int tz_minuteswest; /* minutes W of Greenwich */
    int tz_dsttime; /* type of dst correction */
};

inline int gettimeofday(struct timeval *tv, struct timezone *tz) {
    FILETIME ft;
    unsigned __int64 tmpres = 0;
    static int tzflag;

    if (NULL != tv) {
        GetSystemTimeAsFileTime(&ft);

        tmpres |= ft.dwHighDateTime;
        tmpres <<= 32;
        tmpres |= ft.dwLowDateTime;

        /*converting file time to unix epoch*/
        tmpres /= 10; /*convert into microseconds*/
        tmpres -= DELTA_EPOCH_IN_MICROSECS;
        tv->tv_sec = (long) (tmpres / 1000000UL);
        tv->tv_usec = (long) (tmpres % 1000000UL);
    }

    if (NULL != tz) {
        if (!tzflag) {
            _tzset();
            tzflag++;
        }

        long ltmp;
        _get_timezone(&ltmp);
        tz->tz_minuteswest = ltmp / 60;

        //tz->tz_dsttime = _daylight;
        int itmp;
        _get_daylight(&itmp);
        tz->tz_dsttime = itmp;
    }

    return 0;
}


#endif
