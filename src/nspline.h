/*
 * nspline.h
 *
 * a cubic spline interpolation C library without external dependencies
 *
 * ---------------------------------------------------------------------
 * Copyright (C) 2011, 2014 Tino Kluge (ttk448 at gmail.com)  (original author (C++))
 * Copyright (C) 2019       Turki. s.  (nilputs at gmail.com) (ported to C)
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 *
 */
#ifndef NSPLINE_H
#define NSPLINE_H
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

//define NSP_DEBUG to enable debug assertions
#ifdef NSP_DEBUG
#include <assert.h>
    #define NSP_ASSERT(cond, msg) assert(cond)
#else
    #define NSP_ASSERT(cond, msg) 
#endif

#ifndef NSP_RSCT
    #if !defined(_MSC_VER)
        #define NSP_RSCT restrict
    #endif
    #ifndef NSP_RSCT
        #define NSP_RSCT
    #endif
#endif

#define NMAT_NUPPER 1
#define NMAT_NLOWER 1

struct nspline;
struct nsp_view;

/*--public functions--*/

//returns NSP_OK on success
//preconditions: 
//      * data is sorted by X, so that it's strictly increasing 
//      * there are n >= 3 datapoints
static int nspline_init(struct nspline *ns, struct nsp_view dataview);

//when initialized with this, it doesn't copy, modify nor own the data, but the data must outlive it
static struct nsp_view nsp_const_dview(double *xs, double *ys, long len);
//when initialized with this, a copy is made and the library manages its lifetime
static struct nsp_view nsp_copy_dview(double *xs, double *ys, long len);

//free memory resources held by ns->* 
//doesn't try to free() ns itself, as it doesn't have to be dynamically allocated
static void nspline_deinit(struct nspline *ns);

static double nspline_interpolate(struct nspline *ns, double x);
static double nspline_deriv(struct nspline *ns, int order, double x);

/*--------------------*/

enum NSP_ERR {
    NSP_OK = 0,
    NSP_ALLOC_ERR,
    NSP_INVALID_DVIEW,
    NSP_TOO_FEW, //need at least 3 elements
};

//holds pointers to the dataset
struct nsp_view {
    double *xs;
    double *ys;
    long len;
    bool owns;
};

//a dynamically allocated array of doubles
struct nsp_ddarray {
    double *v;
    long len;
};

//a band matrix with k1 = k2 = 1
struct nmat {
    struct nsp_ddarray bands[3]; //b10 is in bands[0], b00 is in bands[1], b01 is in bands[2].
};

enum NSP_BND_CND {
    NSP_FIRST_DERIV = 1,
    NSP_SECOND_DERIV,
};
enum NSPC {
    NSP_COEFF_A = 0,
    NSP_COEFF_B,
    NSP_COEFF_C,
};

struct nsp_opts {
    bool linearly_extrapolate;
    int left_bnd_cnd;
    int right_bnd_cnd;
    double leftval;
    double rightval;
};

//main datastructure of the library
struct nspline {
    struct nsp_ddarray coeff; //coefficents, stored in a single array as [A, B, C], [A, B, C],
    struct nsp_view dv; //where we store the wild pointers
    struct nsp_opts opts;

    double b0; // interpolation parameters
    double c0; // f(x) = a*(x-x_i)^3 + b*(x-x_i)^2 + c*(x-x_i) + y_i
};

static long nsp_imin(long a, long b) {
    return a < b ? a : b;
}
static long nsp_imax(long a, long b) {
    return a > b ? a : b;
}

static struct nsp_view nsp_const_dview (double *xs, double *ys, long len) {
    struct nsp_view rv = {xs, ys, len, false};
    return rv;
}

//internal detail: on failure returns .xs == .ys == NULL
static struct nsp_view nsp_copy_dview (double *xs, double *ys, long len) {
    struct nsp_view rv = {NULL, NULL, 0, true};
    rv.xs = malloc(sizeof xs[0] * len * 2);
    if (!rv.xs)
        goto fail;
    rv.ys = rv.xs + len;
    memcpy(rv.xs, xs, sizeof xs[0] * len);
    memcpy(rv.ys, ys, sizeof ys[0] * len);
    rv.len = len;
    rv.owns = true;
fail:
    return rv;
}

//internal: this only tells if copying failed
static int nsp_is_dataview_valid__(struct nsp_view view) {
    return view.xs && view.ys && view.len > 0 ? NSP_OK : NSP_INVALID_DVIEW;
}

//internal function, use nsp_const_dview() and nsp_copy_dview() as temporary function arguments
//to avoid double freeing, (nspline_init() calls this, even if it fails)
static void nsp_dataview_deinit(struct nsp_view *view) {
    if (view->owns) {
        NSP_ASSERT(view->ys == view->xs + view->len, "invalid pointer to dataset");
        free(view->xs); //if we allocated them, they're both held by this
    }
    memset(view, 0, sizeof *view);
}

static int nsp_ddarray_init(struct nsp_ddarray *d, long len) {
    d->v = malloc(len * sizeof d->v[0]);
    if (!d->v)
        return NSP_ALLOC_ERR;
    d->len = len;
    return NSP_OK;
}
static void nsp_ddarray_deinit(struct nsp_ddarray *d) {
    free(d->v);
    memset(d, 0, sizeof *d);
}

static int nmat_init(struct nmat *nm, long dim) {
    int rv;
    for (int i=0; i<3; i++) {
        if ((rv = nsp_ddarray_init(&nm->bands[i], dim)) != NSP_OK) {
            for (int j=0; j<i; j++)
                nsp_ddarray_deinit(&nm->bands[i]);
            return rv;
        }
    }
    return NSP_OK;
}

static void nmat_deinit(struct nmat *nm) {
    for (long i=0; i<3; i++) {
        nsp_ddarray_deinit(&nm->bands[i]);
    }
}

static long nmat_dim(struct nmat *nm) {
    NSP_ASSERT(nm->bands[0].len == nm->bands[1].len &&
               nm->bands[1].len == nm->bands[2].len, "internal structures corruption");
    return nm->bands[0].len;
}
//get a pointer to mat[i][j]
static double *nmat_gp(struct nmat *nm, long i, long j) {
    int k = j - i;
    NSP_ASSERT(i >= 0 && i < nmat_dim(nm) && j >= 0 && (k >= -1 && k <= 1), "invalid access in sparse matrix");
    return nm->bands[k+1].v + i;
}
//get mat[i][j]
static double nmat_g(struct nmat *nm, long i, long j) {
    return *nmat_gp(nm, i, j);
}
//set mat[i][j]
static void nmat_s(struct nmat *nm, long i, long j, double value) {
    *nmat_gp(nm, i, j) = value;
}

static struct nsp_opts nsp_default_opts() {
    struct nsp_opts opts = {true, NSP_SECOND_DERIV, NSP_SECOND_DERIV, 0.0, 0.0};
    return opts;
}

//fwd
static int nspline_set_points__(struct nspline *ns);

static int nspline_init(struct nspline *ns, struct nsp_view dataview) {
    int rv = nsp_is_dataview_valid__(dataview);
    if (rv != NSP_OK)
        return rv;
    if (dataview.len < 3) {
        rv = NSP_TOO_FEW;
        goto fail_dataview;
    }
    ns->dv = dataview;
    ns->opts = nsp_default_opts();
    //see struct nspline to see how coefficents are stored
    if ((rv = nsp_ddarray_init(&ns->coeff, dataview.len * 3)) != NSP_OK)
        goto fail_dataview;
    if ((rv = nspline_set_points__(ns)) != NSP_OK)
        goto fail_coeff;

    return NSP_OK;

fail_coeff:
    nsp_ddarray_deinit(&ns->coeff);
fail_dataview:
    nsp_dataview_deinit(&ns->dv);
    return rv;
}

static void nspline_deinit(struct nspline *ns) {
    nsp_ddarray_deinit(&ns->coeff);
    nsp_dataview_deinit(&ns->dv);
}

//get a pointer to a coefficent
static double *nspline_coeff_gp__(struct nspline *ns, enum NSPC which_coeff, long i) {
    NSP_ASSERT(which_coeff >= NSP_COEFF_A && which_coeff <= NSP_COEFF_C, "");
    NSP_ASSERT(i >= 0 && i < ns->dv.len && (ns->coeff.len == ns->dv.len * 3), "");
    return ns->coeff.v + (i * 3) + which_coeff;
}
static double nspline_coeff_g__(struct nspline *ns, enum NSPC which_coeff, long i) {
    return *nspline_coeff_gp__(ns, which_coeff, i);
}
static void nspline_coeff_s__(struct nspline *ns, enum NSPC which_coeff, long i, double value) {
    *nspline_coeff_gp__(ns, which_coeff, i) = value;
}

static void nmat_sol_bbb_0b0_repack(struct nsp_ddarray *in_out, long npacked) {
    NSP_ASSERT(npacked > 0, "");
    long dest = (npacked - 1) * 3 + NSP_COEFF_B;
    NSP_ASSERT(dest < in_out->len, "");
    //[B1B2B3] [B4..] -> [A1B1C1] [A2B2C2] [A3B3C3] [A4B4C4] ...
    for (long i=npacked-1; i>=0; i--) {
        in_out->v[dest] = in_out->v[i];
        dest -= 3;
    }
    NSP_ASSERT(dest + 3 == NSP_COEFF_B, "");
}

static void nmat_lu_decompose(struct nmat *nm, struct nsp_ddarray *lu_array_out) {
    long i_max, j_max, j_min;
    long dim = nmat_dim(nm);
    double x;
    double *lu_arr = lu_array_out->v;

    NSP_ASSERT(lu_array_out->len >= dim, "");
    // preconditioning, normalize column i so that a_ii=1
    for (long i=0; i<dim; i++) {
        NSP_ASSERT(nmat_g(nm, i,i) != 0.0, "");
        lu_arr[i] = 1.0 / nmat_g(nm, i, i);
        j_min = nsp_imax(0, i - NMAT_NLOWER);
        j_max = nsp_imin(dim - 1, i + NMAT_NUPPER);
        for (long j=j_min; j<=j_max; j++) {
            nmat_s(nm, i, j,     nmat_g(nm, i, j) * lu_arr[i]);
        }
        nmat_s(nm, i, i,     1.0); // prevents rounding errors
    }

    // Gauss LR-Decomposition
    for (long k=0; k<dim; k++) {
        i_max = nsp_imin(dim - 1, k + NMAT_NLOWER); // num_lower not a mistake!
        for (long i=k+1; i<=i_max; i++) {
            NSP_ASSERT(nmat_g(nm, k, k) != 0.0, "");
            x = - nmat_g(nm, i, k) / nmat_g(nm, k, k);
            nmat_s(nm, i, k,     -x); //(assembly part of L)
            j_max = nsp_imin(dim - 1, k + NMAT_NUPPER);
            for (long j=k+1; j<=j_max; j++) {
                nmat_s(nm, i, j,     nmat_g(nm, i, j) + x * nmat_g(nm, k, j)); //assembly part of R
            }
        }
    }
}

//internal detail: output must be an already initialized darray with size == nmat_dim(nm)
static int nmat_lu_solve(struct nmat *nm,
                          struct nsp_ddarray *lu_array,
                          struct nsp_ddarray *rhs,
                          struct nsp_ddarray *output)
{
    int rv;
    long dim = nmat_dim(nm);
    NSP_ASSERT(rhs->len == dim, "rhs does not have the correct size");
    NSP_ASSERT(lu_array->len == dim, "lu_array does not have the correct size");
    NSP_ASSERT(output->len >= dim, "output is not initialized with the correct size");

    struct nsp_ddarray tmp;
    if ((rv = nsp_ddarray_init(&tmp, dim)) != NSP_OK) {
        return rv;
    }

    double * NSP_RSCT tmp_arr = tmp.v;
    double * NSP_RSCT rhs_arr = rhs->v;
    double * NSP_RSCT lu_arr  = lu_array->v;
    double * NSP_RSCT out_arr = output->v;

    //y = l_solve(b)
    double sum;
    for (long i=0; i<dim; i++) {
        long j_start = nsp_imax(0, i - NMAT_NLOWER);
        sum = 0;
        for (long j = j_start; j<i; j++)
            sum += nmat_g(nm,i,j) * tmp_arr[j];
        tmp_arr[i] = (rhs_arr[i] * lu_arr[i]) - sum;
    }

    //x = r_solve(y)
    for (long i=dim-1; i>=0; i--) {
        long j_stop = nsp_imin(dim - 1, i + NMAT_NUPPER);
        sum = 0;
        for (long j=i+1; j<=j_stop; j++)
            sum += nmat_g(nm, i, j) * out_arr[j];
        out_arr[i] = (tmp_arr[i] - sum) / nmat_g(nm, i, i);
    }
    
    nsp_ddarray_deinit(&tmp);
    return NSP_OK;
}

//temporary macros that eventually get undef'd
#define GET_COEFF(which, idx)      nspline_coeff_g__(ns, NSP_COEFF_ ## which, idx)
#define SET_COEFF(which, idx, val) nspline_coeff_s__(ns, NSP_COEFF_ ## which, idx, val);

static int nspline_set_points__(struct nspline *ns) {
    NSP_ASSERT(ns->dv.len > 2, "datapoints must be at least 3");
    int rv;
    long n = ns->dv.len;
    double *x = ns->dv.xs;
    double *y = ns->dv.ys;

    struct nmat nm;
    if ((rv = nmat_init(&nm, n)) != NSP_OK)
        return rv;
    struct nsp_ddarray rhs;
    if ((rv = nsp_ddarray_init(&rhs, n)) != NSP_OK)
        goto fail_nmat;
    struct nsp_ddarray lu_array;
    if ((rv = nsp_ddarray_init(&lu_array, n)) != NSP_OK)
        goto fail_rhs;

    #ifdef NSP_DEBUG
        for (long i=0; i < n - 1; i++) {
            NSP_ASSERT(ns->dv.xs[i] < ns->dv.xs[i+1],
                        "data not sorted by x, note: there can't be duplicate x values");
        }
    #endif

    // cubic spline interpolation
    // setting up the matrix and right hand side of the equation system
    // for the parameters b[]
    for (long i=1; i<n-1; i++) {
        nmat_s(&nm, i, i-1,      1.0 / 3.0 * (x[i]   - x[i-1]));
        nmat_s(&nm, i, i,        2.0 / 3.0 * (x[i+1] - x[i-1]));
        nmat_s(&nm, i, i+1,      1.0 / 3.0 * (x[i+1] - x[i]));
        rhs.v[i] = (y[i+1] - y[i]) / (x[i+1]-x[i]) - (y[i] - y[i-1]) / (x[i] - x[i-1]);
    }
    // boundary conditions
    switch (ns->opts.left_bnd_cnd) {
        case NSP_SECOND_DERIV:
            nmat_s(&nm, 0, 0,        2.0);  // 2*b[0] = f''
            nmat_s(&nm, 0, 1,        0.0); 
            rhs.v[0] = ns->opts.leftval;
            break;
        case NSP_FIRST_DERIV:
            // c[0] = f', needs to be re-expressed in terms of b:
            // (2b[0]+b[1])(x[1]-x[0]) = 3 ((y[1]-y[0])/(x[1]-x[0]) - f')
            nmat_s(&nm, 0, 0,        2.0 * (x[1] - x[0]));
            nmat_s(&nm, 0, 1,        1.0 * (x[1] - x[0]));
            rhs.v[0] = 3.0 * ((y[1] - y[0]) / (x[1] - x[0]) - ns->opts.leftval);
            break;
        default:
            NSP_ASSERT(false, "invalid boundary condition");
            break;
    }
    switch (ns->opts.right_bnd_cnd) {
        case NSP_SECOND_DERIV:
            nmat_s(&nm, n-1, n-1,    2.0); // 2*b[n-1] = f''
            nmat_s(&nm, n-1, n-2,    0.0);
            rhs.v[n-1] = ns->opts.rightval;
            break;
        case NSP_FIRST_DERIV:
            // c[n-1] = f', needs to be re-expressed in terms of b: (b[n-2]+2b[n-1])(x[n-1]-x[n-2])
            // = 3 (f' - (y[n-1]-y[n-2])/(x[n-1]-x[n-2]))
            nmat_s(&nm, n-1, n-1,    2.0 * (x[n-1] - x[n-2]));
            nmat_s(&nm, n-1, n-2,    1.0 * (x[n-1] - x[n-2]));
            rhs.v[n-1] = 3.0 * (ns->opts.rightval - (y[n-1] - y[n-2]) / (x[n-1] - x[n-2]));
            break;
        default:
            NSP_ASSERT(false, "invalid boundary condition");
            break;
    }
    nmat_lu_decompose(&nm, &lu_array);
    // solve the equation system to obtain the parameters b[]
    if ((rv = nmat_lu_solve(&nm, &lu_array, &rhs, &ns->coeff)) != NSP_OK)
        goto fail_lu_array;
    //note we passed ns->coeff as the output we now repack it from [B, B, B] to [?,B,?] [?,B,?]
    nmat_sol_bbb_0b0_repack(&ns->coeff, n);
    // calculate parameters a[] and c[] based on b[]
    for (long i=0; i<n-1; i++) {
        double b0 = GET_COEFF(B, i);
        double b1 = GET_COEFF(B, i+1);
        double dx = x[i+1] - x[i];
        double dy = y[i+1] - y[i];
        SET_COEFF(A, i,     1.0 / 3.0 * (b1 - b0) / dx);
        SET_COEFF(C, i,     dy / dx  - 1.0 / 3.0 * (2.0 * b0 + b1) * dx);
    }

    ns->b0 = ns->opts.linearly_extrapolate ? 0.0 : GET_COEFF(B, 0);
    ns->c0 = GET_COEFF(C, 0);

    // for the right extrapolation coefficients
    // f_{n-1}(x) = b*(x-x_{n-1})^2 + c*(x-x_{n-1}) + y_{n-1}
    // m_b[n-1] is determined by the boundary condition
    {
        double am2 = GET_COEFF(A, n-2),  bm2 = GET_COEFF(B, n-2),  cm2 = GET_COEFF(C, n-2);
        double h = x[n-1] - x[n-2];
        SET_COEFF(A, n-1,       0.0);
        SET_COEFF(C, n-1,       3.0 * am2 * h*h + 2.0 * bm2 * h + cm2); // = f'_{n-2}(x_{n-1})
    }
    if (ns->opts.linearly_extrapolate)
        SET_COEFF(B, n-1,       0.0);

    rv = NSP_OK;
fail_lu_array:
    nsp_ddarray_deinit(&lu_array);
fail_rhs:
    nsp_ddarray_deinit(&rhs);
fail_nmat:
    nmat_deinit(&nm);
    return rv;
}

//finds a A[i] <= x, such that A[i+1] is > x
//exclusive end index
static long nsp_bsearch_low_bound(double *arr, long beg, long end, double x) {
    long m = beg + (end - beg) / 2;
    long len = end;
    while (beg < end - 1) {
        if ((arr[m] <= x) && ((m == len - 1) || ((m != len - 1) && arr[m+1] > x))) 
            break;
        else if (arr[m] > x) 
            end = m;
        else
            beg = m;

        m = beg + (end - beg) / 2;
    }
    if (m < end)
        return m;
    return -1;
}

static double nspline_interpolate(struct nspline *ns, double x) {
    long n = ns->dv.len;
    // find the closest point m_x[idx] < x, idx=0 even if x < m_x[0]
    long idx = nsp_bsearch_low_bound(ns->dv.xs, 0, ns->dv.len, x);
    idx = idx == -1 ? 0 : idx;

    double h = x - ns->dv.xs[idx];
    double interpol;
    if (x < ns->dv.xs[0]) { // extrapolation to the left
        interpol = (ns->b0*h + ns->c0) * h + ns->dv.ys[0];
    } 
    else if (x > ns->dv.xs[n-1]) { // extrapolation to the right
        interpol = (GET_COEFF(B, n-1) * h + GET_COEFF(C, n-1)) * h + ns->dv.ys[n-1];
    }
    else { // interpolation
        interpol = ((GET_COEFF(A, idx) * h + GET_COEFF(B, idx)) * h + GET_COEFF(C, idx)) * h + ns->dv.ys[idx];
    }
    return interpol;
}

static double nspline_deriv(struct nspline *ns, int order, double x) {
    NSP_ASSERT(order > 0, "invalid derviative order");
    long n = ns->dv.len;
    long idx = nsp_bsearch_low_bound(ns->dv.xs, 0, ns->dv.len, x);
    idx = idx == -1 ? 0 : idx;

    double h = x - ns->dv.xs[idx];
    double interpol;
    if (x < ns->dv.xs[0]) { // extrapolation to the left
        switch (order) { 
            case 1:  interpol = 2.0 * ns->b0 * h + ns->c0; break;
            case 2:  interpol = 2.0 * ns->b0 * h;          break;
            default: interpol = 0.0;                       break;
        }
    }
    else if (x > ns->dv.xs[n - 1]) { // extrapolation to the right
        switch (order) {
            case 1:  interpol = 2.0 * GET_COEFF(B, n-1) * h + GET_COEFF(C, n-1); break;
            case 2:  interpol = 2.0 * GET_COEFF(B, n-1);                         break;
            default: interpol = 0.0;                                             break;
        }
    }
    else { // interpolation
        double a = GET_COEFF(A, idx),  b = GET_COEFF(B, idx),  c = GET_COEFF(C, idx);
        switch (order) { 
            case 1:  interpol = (3.0 * a * h + 2.0 * b) * h + c; break;
            case 2:  interpol =  6.0 * a * h + 2.0 * b;          break;
            case 3:  interpol =  6.0 * a;                        break;
            default: interpol =  0.0;                            break;
        }
    }
    return interpol;
}

#undef GET_COEFF
#undef SET_COEFF

#endif  //NSPLINE_H
