#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include "nspline.h"
#define MAX_COORD 100

double xs[MAX_COORD];
double ys[MAX_COORD];
int nvals = 0;

void die(const char *msg) {
    fprintf(stderr, "%s\n", msg);
    exit(1);
}
bool checked_read_double(const char *start, const char **end_out, double *out)
{
    const char *end;
    while (*start && !isdigit(*start))
        start++;
    end = start;
    *out = strtod(start, (char **) &end);
    if (end == start)
        return false;
    if (end_out)
        *end_out = end;
    return true;
}
bool checked_read_int(const char *start, int *out) {
    const char *end = start;
    int v;
    v = strtol(start, (char **)&end, 10);
    if (end == start) 
        return false;
    *out = v;
    return true;
}
void die_usage() {
    die("usage: prog [number of points] [start x] [end x]");
}
int main(int argc, const char **argv) {
    if (argc < 4) {
        die_usage();
    }
    int npoints;
    double start_x;
    double end_x;
    if (!checked_read_int(argv[1], &npoints))
        die_usage();
    if (!checked_read_double(argv[2], NULL, &start_x))
        die_usage();
    if (!checked_read_double(argv[3], NULL, &end_x))
        die_usage();
    char buff[256];
    while (fgets(buff, 256, stdin)) {
        const char *start;
        const char *end;
        double x;
        double y;
        start = buff;
        if (!checked_read_double(start, &end, &x))
            break;
        start = end;
        if (!checked_read_double(start, &end, &y))
            break;
        start = end;
        xs[nvals] = x;
        ys[nvals] = y;
        nvals++;
    }
    if (argc > 4 && strcmp(argv[4], "-v") == 0) {
        printf("read %d values\n", nvals);
        for (int i=0; i<nvals; i++) {
            printf("%f %f\n", xs[i], ys[i]);
        }
        printf("-----------------\n");
    }
    struct nspline ns;
    if (nspline_init(&ns, nsp_const_dview(xs, ys, nvals)) != NSP_OK)
        die("nspline failed\n");

    double dx = (end_x - start_x) / npoints;
    while (start_x <= end_x) {
        double value = nspline_interpolate(&ns, start_x);
        printf("%f %f\n", start_x, value);
        start_x += dx;
    }
    nspline_deinit(&ns);

}
