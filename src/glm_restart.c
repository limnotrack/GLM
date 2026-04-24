/******************************************************************************
 *                                                                            *
 * glm_restart.c                                                              *
 *                                                                            *
 * NetCDF restart file writer / reader for GLM.                               *
 *                                                                            *
 * State saved:                                                               *
 *   - Lake layer arrays (Temp, Salinity, Height, LayerVol, LayerArea,        *
 *     MeanHeight, Density, Epsilon, Umean, Uorb, LayerStress)               *
 *   - WQ per-layer variables (WQ_Vars)                                       *
 *   - Surface state (AvgSurfTemp, delzBlueIce, delzWhiteIce, delzSnow,      *
 *     RhoSnow)                                                               *
 *   - Mixer state (DepMX, PrevThick, Energy_AvailableMix, OldSlope,         *
 *     Time_count_end_shear, Time_count_sim, Half_Seiche_Period, FO, FSUM,   *
 *     u_f, u0, u_avg, Mixer_Count)                                           *
 *   - Benthic WQ sheet variables (WQS_Vars, Num_WQ_Ben)                     *
 *   - Per-inflow insertion queue (DOld, QIns, TIns, SIns, DIIns, InPar,     *
 *     NoIns, iCnt, SubmElev, WQIns)                                         *
 *   - Sediment layer temps per zone when sed_heat_model == 2                *
 *   - Particle (PTM) state: PTM_Stat and PTM_Vars (when ptm_sw is TRUE)     *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Agriculture and Environment                                  *
 *     The University of Western Australia                                    *
 *                                                                            *
 * Copyright 2013-2026 : The University of Western Australia                  *
 *                                                                            *
 *  This file is part of GLM (General Lake Model)                             *
 *                                                                            *
 *  GLM is free software: you can redistribute it and/or modify               *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  GLM is distributed in the hope that it will be useful,                    *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
 *                                                                            *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <netcdf.h>

#include "glm.h"
#include "glm_types.h"
#include "glm_globals.h"
#include "glm_ncdf.h"
#include "glm_restart.h"

/*----------------------------------------------------------------------------*/
/* Global restart configuration                                               */
char *rst_fn    = NULL;
int   rst_nsave = 0;

/*----------------------------------------------------------------------------*/
/* Helpers                                                                    */
/*----------------------------------------------------------------------------*/

#define RST_CHECK(call) do { \
    int _err = (call); \
    if (_err != NC_NOERR) { \
        fprintf(stderr, "NetCDF restart error at %s:%d : %s\n", \
                __FILE__, __LINE__, nc_strerror(_err)); \
        nc_close(rst_ncid_); \
        exit(1); \
    } \
} while (0)

/* Helper that aborts with an error including context. Used only during open
   so we can close safely.  We keep a file-scope variable set before each
   macro block so the RST_CHECK macro can close on error. */
static int rst_ncid_ = -1;

static void def_var_d(int ncid, const char *name, int ndim, const int *dims,
                      int *id)
{
    RST_CHECK(nc_def_var(ncid, name, NC_REALTYPE, ndim, dims, id));
}

static void def_var_i(int ncid, const char *name, int ndim, const int *dims,
                      int *id)
{
    RST_CHECK(nc_def_var(ncid, name, NC_INT, ndim, dims, id));
}

/*----------------------------------------------------------------------------*/
/* write_glm_restart                                                          */
/*----------------------------------------------------------------------------*/
void write_glm_restart(const char *fn)
{
    int ncid;
    rst_ncid_ = -1;

    /* ---- create/overwrite the file ---- */
    RST_CHECK(nc_create(fn, NC_CLOBBER|NC_NETCDF4, &ncid));
    rst_ncid_ = ncid;

    /* -------- define dimensions -------- */
    int dim_nlev, dim_maxpar, dim_numinf, dim_wqvars, dim_wqben, dim_zones;
    int dim_sedlayers = -1;

    RST_CHECK(nc_def_dim(ncid, "nlev",        NumLayers,   &dim_nlev));
    RST_CHECK(nc_def_dim(ncid, "max_par",     MaxPar,      &dim_maxpar));
    RST_CHECK(nc_def_dim(ncid, "num_inflows",  NumInf > 0 ? NumInf : 1,
                         &dim_numinf));
    int nwq = (wq_calc && Tot_WQ_Vars > 0) ? Tot_WQ_Vars : 1;
    RST_CHECK(nc_def_dim(ncid, "wq_vars",     nwq,         &dim_wqvars));
    int nben = (wq_calc && Num_WQ_Ben > 0) ? Num_WQ_Ben : 1;
    RST_CHECK(nc_def_dim(ncid, "wq_ben",      nben,        &dim_wqben));
    RST_CHECK(nc_def_dim(ncid, "n_zones",     n_zones > 0 ? n_zones : 1,
                         &dim_zones));
    int n_sed_layers_rst = 0;
    if (sed_heat_model == 2 && n_zones > 0 && theZones != NULL)
        n_sed_layers_rst = theZones[0].n_sed_layers;
    if (n_sed_layers_rst > 0)
        RST_CHECK(nc_def_dim(ncid, "sed_layers", n_sed_layers_rst,
                             &dim_sedlayers));

    /* PTM dimensions (only defined when particles are active) */
    static const int PTM_STAT_NVARS = 6; /* STAT,IDX2,IDX3,LAYR,FLAG,PTID */
    int dim_ptm_par = -1, dim_ptm_sv = -1, dim_ptm_wqv = -1;
    int ptm_enabled = (ptm_sw && PTM_Stat != NULL && max_particle_num > 0) ? 1 : 0;
    if (ptm_enabled) {
        RST_CHECK(nc_def_dim(ncid, "ptm_particles", max_particle_num, &dim_ptm_par));
        RST_CHECK(nc_def_dim(ncid, "ptm_stat_vars", PTM_STAT_NVARS,   &dim_ptm_sv));
        RST_CHECK(nc_def_dim(ncid, "ptm_wq_vars",   Num_WQ_Vars,      &dim_ptm_wqv));
    }

    /* -------- global attributes -------- */
    RST_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "NumLayers",       NC_INT, 1, &NumLayers));
    RST_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "NumInf",          NC_INT, 1, &NumInf));
    RST_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "Tot_WQ_Vars",     NC_INT, 1, &Tot_WQ_Vars));
    RST_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "Num_WQ_Ben",      NC_INT, 1, &Num_WQ_Ben));
    RST_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "n_zones",         NC_INT, 1, &n_zones));
    RST_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "sed_heat_model",  NC_INT, 1, &sed_heat_model));
    RST_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "n_sed_layers_rst",NC_INT, 1, &n_sed_layers_rst));
    RST_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "ptm_enabled",     NC_INT, 1, &ptm_enabled));
    RST_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "max_particle_num",NC_INT, 1, &max_particle_num));
    RST_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "Num_WQ_Vars",     NC_INT, 1, &Num_WQ_Vars));

    /* -------- define variables -------- */

    /* Lake layer arrays [nlev] */
    int id_temp, id_salt, id_height, id_layvol, id_layarea;
    int id_meanht, id_dens, id_eps, id_umean, id_uorb, id_stress;
    def_var_d(ncid, "lake_temp",        1, &dim_nlev, &id_temp);
    def_var_d(ncid, "lake_salt",        1, &dim_nlev, &id_salt);
    def_var_d(ncid, "lake_height",      1, &dim_nlev, &id_height);
    def_var_d(ncid, "lake_layervol",    1, &dim_nlev, &id_layvol);
    def_var_d(ncid, "lake_layerarea",   1, &dim_nlev, &id_layarea);
    def_var_d(ncid, "lake_meanheight",  1, &dim_nlev, &id_meanht);
    def_var_d(ncid, "lake_density",     1, &dim_nlev, &id_dens);
    def_var_d(ncid, "lake_epsilon",     1, &dim_nlev, &id_eps);
    def_var_d(ncid, "lake_umean",       1, &dim_nlev, &id_umean);
    def_var_d(ncid, "lake_uorb",        1, &dim_nlev, &id_uorb);
    def_var_d(ncid, "lake_stress",      1, &dim_nlev, &id_stress);

    /* WQ per-layer [wq_vars, nlev] (stored as flat [wq_vars * nlev]) */
    int id_wq = -1;
    if (wq_calc && Tot_WQ_Vars > 0) {
        int dims2[2] = { dim_wqvars, dim_nlev };
        def_var_d(ncid, "wq_vars", 2, dims2, &id_wq);
    }

    /* Benthic WQ sheet [wq_ben] */
    int id_wqs = -1;
    if (wq_calc && Num_WQ_Ben > 0 && WQS_Vars != NULL) {
        def_var_d(ncid, "wqs_vars", 1, &dim_wqben, &id_wqs);
    }

    /* Surface scalars */
    int id_avgt, id_blue, id_white, id_snow, id_rhosnow;
    def_var_d(ncid, "avg_surf_temp",    0, NULL, &id_avgt);
    def_var_d(ncid, "blue_ice",         0, NULL, &id_blue);
    def_var_d(ncid, "white_ice",        0, NULL, &id_white);
    def_var_d(ncid, "snow_thickness",   0, NULL, &id_snow);
    def_var_d(ncid, "rho_snow",         0, NULL, &id_rhosnow);

    /* Mixer scalars */
    int id_depmx, id_prevthick, id_eavail, id_oldslope;
    int id_tce, id_tcs, id_hsp, id_fo, id_fsum, id_uf, id_u0, id_uavg;
    int id_mxcnt;
    def_var_d(ncid, "mixer_dep_mx",            0, NULL, &id_depmx);
    def_var_d(ncid, "mixer_prev_thick",        0, NULL, &id_prevthick);
    def_var_d(ncid, "mixer_energy_available",  0, NULL, &id_eavail);
    def_var_d(ncid, "mixer_old_slope",         0, NULL, &id_oldslope);
    def_var_d(ncid, "mixer_time_count_end_shear", 0, NULL, &id_tce);
    def_var_d(ncid, "mixer_time_count_sim",    0, NULL, &id_tcs);
    def_var_d(ncid, "mixer_half_seiche_period",0, NULL, &id_hsp);
    def_var_d(ncid, "mixer_fo",                0, NULL, &id_fo);
    def_var_d(ncid, "mixer_fsum",              0, NULL, &id_fsum);
    def_var_d(ncid, "mixer_u_f",               0, NULL, &id_uf);
    def_var_d(ncid, "mixer_u0",                0, NULL, &id_u0);
    def_var_d(ncid, "mixer_u_avg",             0, NULL, &id_uavg);
    def_var_i(ncid, "mixer_count",             0, NULL, &id_mxcnt);

    /* Per-inflow insertion queue [num_inflows, max_par] */
    int id_dold, id_qins, id_tins, id_sins, id_diins;
    int id_inpar, id_noins, id_icnt, id_submelev;
    int id_wqins = -1;
    if (NumInf > 0) {
        int dims_ip[2] = { dim_numinf, dim_maxpar };
        def_var_d(ncid, "inflow_dold",  2, dims_ip, &id_dold);
        def_var_d(ncid, "inflow_qins",  2, dims_ip, &id_qins);
        def_var_d(ncid, "inflow_tins",  2, dims_ip, &id_tins);
        def_var_d(ncid, "inflow_sins",  2, dims_ip, &id_sins);
        def_var_d(ncid, "inflow_diins", 2, dims_ip, &id_diins);
        def_var_i(ncid, "inflow_inpar", 2, dims_ip, &id_inpar);
        def_var_i(ncid, "inflow_noins", 1, &dim_numinf, &id_noins);
        def_var_i(ncid, "inflow_icnt",  1, &dim_numinf, &id_icnt);
        def_var_d(ncid, "inflow_submelev", 1, &dim_numinf, &id_submelev);
        if (wq_calc && Tot_WQ_Vars > 0) {
            int dims_iw[3] = { dim_numinf, dim_maxpar, dim_wqvars };
            def_var_d(ncid, "inflow_wqins", 3, dims_iw, &id_wqins);
        }
    }

    /* Sediment layer temps per zone [n_zones, sed_layers] when model==2 */
    int id_sedtemp = -1, id_sedhf = -1;
    if (n_sed_layers_rst > 0) {
        int dims_sz[2] = { dim_zones, dim_sedlayers };
        def_var_d(ncid, "sed_layer_temp", 2, dims_sz, &id_sedtemp);
        def_var_d(ncid, "sed_heatflux",   1, &dim_zones, &id_sedhf);
    }

    /* PTM particle state [ptm_stat_vars, ptm_particles] and
     *                    [ptm_wq_vars,   ptm_particles]             */
    int id_ptm_stat = -1, id_ptm_vars = -1;
    if (ptm_enabled) {
        int dims_ps[2] = { dim_ptm_sv,  dim_ptm_par };
        int dims_pv[2] = { dim_ptm_wqv, dim_ptm_par };
        def_var_i(ncid, "ptm_stat", 2, dims_ps, &id_ptm_stat);
        def_var_d(ncid, "ptm_vars", 2, dims_pv, &id_ptm_vars);
    }

    /* End define mode */
    RST_CHECK(nc_enddef(ncid));

    /* -------- write data -------- */

    /* Allocate a working buffer large enough for any 1-D layer array */
    AED_REAL *buf = malloc(NumLayers * sizeof(AED_REAL));
    if (!buf) { fprintf(stderr, "glm_restart: out of memory\n"); exit(1); }

    /* Lake layer fields */
#define WRITE_LAKE_FIELD(varid, field) \
    do { for (int _i = 0; _i < NumLayers; _i++) buf[_i] = Lake[_i].field; \
         RST_CHECK(nc_put_var_double(ncid, varid, buf)); } while(0)

    WRITE_LAKE_FIELD(id_temp,    Temp);
    WRITE_LAKE_FIELD(id_salt,    Salinity);
    WRITE_LAKE_FIELD(id_height,  Height);
    WRITE_LAKE_FIELD(id_layvol,  LayerVol);
    WRITE_LAKE_FIELD(id_layarea, LayerArea);
    WRITE_LAKE_FIELD(id_meanht,  MeanHeight);
    WRITE_LAKE_FIELD(id_dens,    Density);
    WRITE_LAKE_FIELD(id_eps,     Epsilon);
    WRITE_LAKE_FIELD(id_umean,   Umean);
    WRITE_LAKE_FIELD(id_uorb,    Uorb);
    WRITE_LAKE_FIELD(id_stress,  LayerStress);
#undef WRITE_LAKE_FIELD

    free(buf);

    /* WQ per-layer: WQ_Vars layout is [Tot_WQ_Vars * MaxLayers] with index
     * _WQ_Vars(var, lyr) = WQ_Vars[Tot_WQ_Vars * lyr + var]
     * We write a [wq_vars, nlev] array (contiguous in var for each layer). */
    if (id_wq >= 0) {
        /* Build a compact [Tot_WQ_Vars * NumLayers] buffer */
        AED_REAL *wqbuf = malloc(Tot_WQ_Vars * NumLayers * sizeof(AED_REAL));
        if (!wqbuf) { fprintf(stderr, "glm_restart: out of memory\n"); exit(1); }
        for (int v = 0; v < Tot_WQ_Vars; v++)
            for (int l = 0; l < NumLayers; l++)
                wqbuf[v * NumLayers + l] = _WQ_Vars(v, l);
        RST_CHECK(nc_put_var_double(ncid, id_wq, wqbuf));
        free(wqbuf);
    }

    /* Benthic WQ sheet */
    if (id_wqs >= 0)
        RST_CHECK(nc_put_var_double(ncid, id_wqs, WQS_Vars));

    /* Surface scalars */
    RST_CHECK(nc_put_var_double(ncid, id_avgt,    &AvgSurfTemp));
    RST_CHECK(nc_put_var_double(ncid, id_blue,    &SurfData.delzBlueIce));
    RST_CHECK(nc_put_var_double(ncid, id_white,   &SurfData.delzWhiteIce));
    RST_CHECK(nc_put_var_double(ncid, id_snow,    &SurfData.delzSnow));
    RST_CHECK(nc_put_var_double(ncid, id_rhosnow, &SurfData.RhoSnow));

    /* Mixer scalars */
    RST_CHECK(nc_put_var_double(ncid, id_depmx,    &DepMX));
    RST_CHECK(nc_put_var_double(ncid, id_prevthick, &PrevThick));
    RST_CHECK(nc_put_var_double(ncid, id_eavail,   &Energy_AvailableMix));
    RST_CHECK(nc_put_var_double(ncid, id_oldslope, &OldSlope));
    RST_CHECK(nc_put_var_double(ncid, id_tce,      &Time_count_end_shear));
    RST_CHECK(nc_put_var_double(ncid, id_tcs,      &Time_count_sim));
    RST_CHECK(nc_put_var_double(ncid, id_hsp,      &Half_Seiche_Period));
    RST_CHECK(nc_put_var_double(ncid, id_fo,       &FO));
    RST_CHECK(nc_put_var_double(ncid, id_fsum,     &FSUM));
    RST_CHECK(nc_put_var_double(ncid, id_uf,       &u_f));
    RST_CHECK(nc_put_var_double(ncid, id_u0,       &u0));
    RST_CHECK(nc_put_var_double(ncid, id_uavg,     &u_avg));
    RST_CHECK(nc_put_var_int   (ncid, id_mxcnt,    &Mixer_Count));

    /* Per-inflow insertion queue */
    if (NumInf > 0) {
        /* Build flat [NumInf * MaxPar] buffers */
        AED_REAL *dbuf = malloc(NumInf * MaxPar * sizeof(AED_REAL));
        int      *ibuf = malloc(NumInf * MaxPar * sizeof(int));
        AED_REAL *svec = malloc(NumInf         * sizeof(AED_REAL));
        int      *ivec = malloc(NumInf         * sizeof(int));
        if (!dbuf || !ibuf || !svec || !ivec) {
            fprintf(stderr, "glm_restart: out of memory\n"); exit(1);
        }

#define FILL_INF_DFIELD(arr) \
        do { for (int _s = 0; _s < NumInf; _s++) \
                 for (int _p = 0; _p < MaxPar; _p++) \
                     dbuf[_s * MaxPar + _p] = Inflows[_s].arr[_p]; } while(0)

#define FILL_INF_IFIELD(arr) \
        do { for (int _s = 0; _s < NumInf; _s++) \
                 for (int _p = 0; _p < MaxPar; _p++) \
                     ibuf[_s * MaxPar + _p] = Inflows[_s].arr[_p]; } while(0)

        FILL_INF_DFIELD(DOld);
        RST_CHECK(nc_put_var_double(ncid, id_dold, dbuf));
        FILL_INF_DFIELD(QIns);
        RST_CHECK(nc_put_var_double(ncid, id_qins, dbuf));
        FILL_INF_DFIELD(TIns);
        RST_CHECK(nc_put_var_double(ncid, id_tins, dbuf));
        FILL_INF_DFIELD(SIns);
        RST_CHECK(nc_put_var_double(ncid, id_sins, dbuf));
        FILL_INF_DFIELD(DIIns);
        RST_CHECK(nc_put_var_double(ncid, id_diins, dbuf));
        FILL_INF_IFIELD(InPar);
        RST_CHECK(nc_put_var_int(ncid, id_inpar, ibuf));
#undef FILL_INF_DFIELD
#undef FILL_INF_IFIELD

        for (int s = 0; s < NumInf; s++) ivec[s] = Inflows[s].NoIns;
        RST_CHECK(nc_put_var_int(ncid, id_noins, ivec));
        for (int s = 0; s < NumInf; s++) ivec[s] = Inflows[s].iCnt;
        RST_CHECK(nc_put_var_int(ncid, id_icnt, ivec));
        for (int s = 0; s < NumInf; s++) svec[s] = Inflows[s].SubmElev;
        RST_CHECK(nc_put_var_double(ncid, id_submelev, svec));

        /* WQIns [NumInf, MaxPar, Tot_WQ_Vars] */
        if (id_wqins >= 0) {
            AED_REAL *wibuf = malloc(NumInf * MaxPar * Tot_WQ_Vars * sizeof(AED_REAL));
            if (!wibuf) { fprintf(stderr, "glm_restart: out of memory\n"); exit(1); }
            for (int s = 0; s < NumInf; s++)
                for (int p = 0; p < MaxPar; p++)
                    for (int v = 0; v < Tot_WQ_Vars; v++)
                        wibuf[(s * MaxPar + p) * Tot_WQ_Vars + v] =
                            Inflows[s].WQIns[p][v];
            RST_CHECK(nc_put_var_double(ncid, id_wqins, wibuf));
            free(wibuf);
        }

        free(dbuf); free(ibuf); free(svec); free(ivec);
    }

    /* Sediment layer data when sed_heat_model == 2 */
    if (n_sed_layers_rst > 0) {
        AED_REAL *stbuf = malloc(n_zones * n_sed_layers_rst * sizeof(AED_REAL));
        AED_REAL *hfbuf = malloc(n_zones * sizeof(AED_REAL));
        if (!stbuf || !hfbuf) { fprintf(stderr, "glm_restart: out of memory\n"); exit(1); }
        for (int z = 0; z < n_zones; z++) {
            hfbuf[z] = theZones[z].heatflux;
            for (int k = 0; k < n_sed_layers_rst; k++)
                stbuf[z * n_sed_layers_rst + k] = theZones[z].layers[k].temp;
        }
        RST_CHECK(nc_put_var_double(ncid, id_sedtemp, stbuf));
        RST_CHECK(nc_put_var_double(ncid, id_sedhf,   hfbuf));
        free(stbuf); free(hfbuf);
    }

    /* PTM particle state
     * PTM_Stat layout: PTM_Stat[var * max_particle_num + part]  (grp=0 only)
     * PTM_Vars layout: PTM_Vars[var * max_particle_num + part]  (grp=0 only)
     * These match the [ptm_stat_vars/ptm_wq_vars, ptm_particles] NC layout. */
    if (ptm_enabled) {
        RST_CHECK(nc_put_var_int   (ncid, id_ptm_stat, PTM_Stat));
        RST_CHECK(nc_put_var_double(ncid, id_ptm_vars, PTM_Vars));
    }

    nc_close(ncid);
    rst_ncid_ = -1;
}

/*----------------------------------------------------------------------------*/
/* read_glm_restart                                                           */
/*----------------------------------------------------------------------------*/
int read_glm_restart(const char *fn)
{
    int ncid;
    int err;

    err = nc_open(fn, NC_NOWRITE, &ncid);
    if (err == NC_ENOTFOUND) return 0; /* file does not exist */
    if (err != NC_NOERR) {
        /* Extra fallback: if it's some other open error treat absent file gracefully */
        FILE *_f = fopen(fn, "r");
        if (_f == NULL) return 0;
        fclose(_f);
        fprintf(stderr, "glm_restart: cannot open %s : %s\n", fn, nc_strerror(err));
        exit(1);
    }
    rst_ncid_ = ncid;

    /* ---- validate global attributes ---- */
    int att_nlev, att_numinf, att_tot_wq, att_nben, att_nzones;
    int att_sed_heat_model, att_n_sed_layers_rst;
    RST_CHECK(nc_get_att_int(ncid, NC_GLOBAL, "NumLayers",      &att_nlev));
    RST_CHECK(nc_get_att_int(ncid, NC_GLOBAL, "NumInf",         &att_numinf));
    RST_CHECK(nc_get_att_int(ncid, NC_GLOBAL, "Tot_WQ_Vars",    &att_tot_wq));
    RST_CHECK(nc_get_att_int(ncid, NC_GLOBAL, "Num_WQ_Ben",     &att_nben));
    RST_CHECK(nc_get_att_int(ncid, NC_GLOBAL, "n_zones",        &att_nzones));
    RST_CHECK(nc_get_att_int(ncid, NC_GLOBAL, "sed_heat_model", &att_sed_heat_model));
    RST_CHECK(nc_get_att_int(ncid, NC_GLOBAL, "n_sed_layers_rst", &att_n_sed_layers_rst));

    if (att_nlev > MaxLayers) {
        fprintf(stderr, "glm_restart: restart NumLayers (%d) > MaxLayers (%d)\n",
                att_nlev, MaxLayers);
        exit(1);
    }
    if (att_numinf != NumInf) {
        fprintf(stderr, "glm_restart: restart NumInf (%d) != current NumInf (%d)\n",
                att_numinf, NumInf);
        exit(1);
    }
    if (att_tot_wq != Tot_WQ_Vars) {
        fprintf(stderr, "glm_restart: restart Tot_WQ_Vars (%d) != current (%d)\n",
                att_tot_wq, Tot_WQ_Vars);
        exit(1);
    }
    if (att_nben != Num_WQ_Ben) {
        fprintf(stderr, "glm_restart: restart Num_WQ_Ben (%d) != current (%d)\n",
                att_nben, Num_WQ_Ben);
        exit(1);
    }

    NumLayers = att_nlev;

    /* ---- helper lambda-like macro ---- */
#define RD_VAR(name, varid) \
    int varid; \
    RST_CHECK(nc_inq_varid(ncid, name, &varid));

    /* ---- lake layer fields ---- */
    {
        AED_REAL *buf = malloc(NumLayers * sizeof(AED_REAL));
        if (!buf) { fprintf(stderr, "glm_restart: out of memory\n"); exit(1); }

#define READ_LAKE_FIELD(varname, field) \
        do { int _id; RST_CHECK(nc_inq_varid(ncid, varname, &_id)); \
             RST_CHECK(nc_get_var_double(ncid, _id, buf)); \
             for (int _i = 0; _i < NumLayers; _i++) Lake[_i].field = buf[_i]; \
        } while(0)

        READ_LAKE_FIELD("lake_temp",       Temp);
        READ_LAKE_FIELD("lake_salt",       Salinity);
        READ_LAKE_FIELD("lake_height",     Height);
        READ_LAKE_FIELD("lake_layervol",   LayerVol);
        READ_LAKE_FIELD("lake_layerarea",  LayerArea);
        READ_LAKE_FIELD("lake_meanheight", MeanHeight);
        READ_LAKE_FIELD("lake_density",    Density);
        READ_LAKE_FIELD("lake_epsilon",    Epsilon);
        READ_LAKE_FIELD("lake_umean",      Umean);
        READ_LAKE_FIELD("lake_uorb",       Uorb);
        READ_LAKE_FIELD("lake_stress",     LayerStress);
#undef READ_LAKE_FIELD

        free(buf);
    }

    /* ---- WQ per-layer ---- */
    if (wq_calc && Tot_WQ_Vars > 0 && WQ_Vars != NULL) {
        int id_wq;
        RST_CHECK(nc_inq_varid(ncid, "wq_vars", &id_wq));
        AED_REAL *wqbuf = malloc(Tot_WQ_Vars * NumLayers * sizeof(AED_REAL));
        if (!wqbuf) { fprintf(stderr, "glm_restart: out of memory\n"); exit(1); }
        RST_CHECK(nc_get_var_double(ncid, id_wq, wqbuf));
        for (int v = 0; v < Tot_WQ_Vars; v++)
            for (int l = 0; l < NumLayers; l++)
                _WQ_Vars(v, l) = wqbuf[v * NumLayers + l];
        free(wqbuf);
    }

    /* ---- Benthic WQ sheet ---- */
    if (wq_calc && Num_WQ_Ben > 0 && WQS_Vars != NULL) {
        int id_wqs;
        RST_CHECK(nc_inq_varid(ncid, "wqs_vars", &id_wqs));
        RST_CHECK(nc_get_var_double(ncid, id_wqs, WQS_Vars));
    }

    /* ---- Surface scalars ---- */
    {
#define RD_SCALAR_D(name, dst) \
        do { int _id; RST_CHECK(nc_inq_varid(ncid, name, &_id)); \
             RST_CHECK(nc_get_var_double(ncid, _id, &(dst))); } while(0)
#define RD_SCALAR_I(name, dst) \
        do { int _id; RST_CHECK(nc_inq_varid(ncid, name, &_id)); \
             RST_CHECK(nc_get_var_int(ncid, _id, &(dst))); } while(0)

        RD_SCALAR_D("avg_surf_temp",  AvgSurfTemp);
        RD_SCALAR_D("blue_ice",       SurfData.delzBlueIce);
        RD_SCALAR_D("white_ice",      SurfData.delzWhiteIce);
        RD_SCALAR_D("snow_thickness", SurfData.delzSnow);
        RD_SCALAR_D("rho_snow",       SurfData.RhoSnow);

        /* ---- Mixer scalars ---- */
        RD_SCALAR_D("mixer_dep_mx",             DepMX);
        RD_SCALAR_D("mixer_prev_thick",         PrevThick);
        RD_SCALAR_D("mixer_energy_available",   Energy_AvailableMix);
        RD_SCALAR_D("mixer_old_slope",          OldSlope);
        RD_SCALAR_D("mixer_time_count_end_shear", Time_count_end_shear);
        RD_SCALAR_D("mixer_time_count_sim",     Time_count_sim);
        RD_SCALAR_D("mixer_half_seiche_period", Half_Seiche_Period);
        RD_SCALAR_D("mixer_fo",                 FO);
        RD_SCALAR_D("mixer_fsum",               FSUM);
        RD_SCALAR_D("mixer_u_f",                u_f);
        RD_SCALAR_D("mixer_u0",                 u0);
        RD_SCALAR_D("mixer_u_avg",              u_avg);
        RD_SCALAR_I("mixer_count",              Mixer_Count);
#undef RD_SCALAR_D
#undef RD_SCALAR_I
    }

    /* Update ice flag */
    if (SurfData.delzBlueIce > 0.0 || SurfData.delzWhiteIce > 0.0) {
        ice = TRUE;
    }

    /* ---- Per-inflow insertion queue ---- */
    if (NumInf > 0) {
        AED_REAL *dbuf = malloc(NumInf * MaxPar * sizeof(AED_REAL));
        int      *ibuf = malloc(NumInf * MaxPar * sizeof(int));
        AED_REAL *svec = malloc(NumInf         * sizeof(AED_REAL));
        int      *ivec = malloc(NumInf         * sizeof(int));
        if (!dbuf || !ibuf || !svec || !ivec) {
            fprintf(stderr, "glm_restart: out of memory\n"); exit(1);
        }

#define READ_INF_DFIELD(varname, arr) \
        do { int _id; RST_CHECK(nc_inq_varid(ncid, varname, &_id)); \
             RST_CHECK(nc_get_var_double(ncid, _id, dbuf)); \
             for (int _s = 0; _s < NumInf; _s++) \
                 for (int _p = 0; _p < MaxPar; _p++) \
                     Inflows[_s].arr[_p] = dbuf[_s * MaxPar + _p]; } while(0)

#define READ_INF_IFIELD(varname, arr) \
        do { int _id; RST_CHECK(nc_inq_varid(ncid, varname, &_id)); \
             RST_CHECK(nc_get_var_int(ncid, _id, ibuf)); \
             for (int _s = 0; _s < NumInf; _s++) \
                 for (int _p = 0; _p < MaxPar; _p++) \
                     Inflows[_s].arr[_p] = ibuf[_s * MaxPar + _p]; } while(0)

        READ_INF_DFIELD("inflow_dold",  DOld);
        READ_INF_DFIELD("inflow_qins",  QIns);
        READ_INF_DFIELD("inflow_tins",  TIns);
        READ_INF_DFIELD("inflow_sins",  SIns);
        READ_INF_DFIELD("inflow_diins", DIIns);
        READ_INF_IFIELD("inflow_inpar", InPar);
#undef READ_INF_DFIELD
#undef READ_INF_IFIELD

        {
            int _id;
            RST_CHECK(nc_inq_varid(ncid, "inflow_noins", &_id));
            RST_CHECK(nc_get_var_int(ncid, _id, ivec));
            for (int s = 0; s < NumInf; s++) Inflows[s].NoIns = ivec[s];

            RST_CHECK(nc_inq_varid(ncid, "inflow_icnt", &_id));
            RST_CHECK(nc_get_var_int(ncid, _id, ivec));
            for (int s = 0; s < NumInf; s++) Inflows[s].iCnt = ivec[s];

            RST_CHECK(nc_inq_varid(ncid, "inflow_submelev", &_id));
            RST_CHECK(nc_get_var_double(ncid, _id, svec));
            for (int s = 0; s < NumInf; s++) Inflows[s].SubmElev = svec[s];
        }

        if (wq_calc && Tot_WQ_Vars > 0) {
            AED_REAL *wibuf = malloc(NumInf * MaxPar * Tot_WQ_Vars * sizeof(AED_REAL));
            if (!wibuf) { fprintf(stderr, "glm_restart: out of memory\n"); exit(1); }
            int _id;
            RST_CHECK(nc_inq_varid(ncid, "inflow_wqins", &_id));
            RST_CHECK(nc_get_var_double(ncid, _id, wibuf));
            for (int s = 0; s < NumInf; s++)
                for (int p = 0; p < MaxPar; p++)
                    for (int v = 0; v < Tot_WQ_Vars; v++)
                        Inflows[s].WQIns[p][v] =
                            wibuf[(s * MaxPar + p) * Tot_WQ_Vars + v];
            free(wibuf);
        }

        free(dbuf); free(ibuf); free(svec); free(ivec);
    }

    /* ---- Sediment layer data ---- */
    if (att_n_sed_layers_rst > 0 && sed_heat_model == 2 &&
        n_zones > 0 && theZones != NULL) {
        int nsrl = att_n_sed_layers_rst;
        AED_REAL *stbuf = malloc(n_zones * nsrl * sizeof(AED_REAL));
        AED_REAL *hfbuf = malloc(n_zones * sizeof(AED_REAL));
        if (!stbuf || !hfbuf) { fprintf(stderr, "glm_restart: out of memory\n"); exit(1); }
        int _id;
        RST_CHECK(nc_inq_varid(ncid, "sed_layer_temp", &_id));
        RST_CHECK(nc_get_var_double(ncid, _id, stbuf));
        RST_CHECK(nc_inq_varid(ncid, "sed_heatflux", &_id));
        RST_CHECK(nc_get_var_double(ncid, _id, hfbuf));
        for (int z = 0; z < n_zones; z++) {
            theZones[z].heatflux = hfbuf[z];
            int kmax = (nsrl < theZones[z].n_sed_layers) ?
                        nsrl : theZones[z].n_sed_layers;
            for (int k = 0; k < kmax; k++)
                theZones[z].layers[k].temp = stbuf[z * nsrl + k];
        }
        free(stbuf); free(hfbuf);
    }

    /* ---- PTM particle state ---- */
    {
        int att_ptm_en = 0, att_max_ptm = 0, att_nwq = 0;
        /* Use nc_get_att_int without RST_CHECK — attribute may be absent in
         * older restart files; errors here are non-fatal. */
        nc_get_att_int(ncid, NC_GLOBAL, "ptm_enabled",     &att_ptm_en);
        nc_get_att_int(ncid, NC_GLOBAL, "max_particle_num",&att_max_ptm);
        nc_get_att_int(ncid, NC_GLOBAL, "Num_WQ_Vars",     &att_nwq);

        if (att_ptm_en) {
            if (!ptm_sw || PTM_Stat == NULL) {
                fprintf(stderr, "     WARNING: restart has PTM state but "
                        "ptm_sw is off; particle state skipped.\n");
            } else if (att_max_ptm != max_particle_num) {
                fprintf(stderr, "     WARNING: restart max_particle_num (%d) "
                        "!= current (%d); particle state skipped.\n",
                        att_max_ptm, max_particle_num);
            } else if (att_nwq != Num_WQ_Vars) {
                fprintf(stderr, "     WARNING: restart Num_WQ_Vars (%d) "
                        "!= current (%d); particle state skipped.\n",
                        att_nwq, Num_WQ_Vars);
            } else {
                int _id;
                RST_CHECK(nc_inq_varid(ncid, "ptm_stat", &_id));
                RST_CHECK(nc_get_var_int(ncid, _id, PTM_Stat));
                RST_CHECK(nc_inq_varid(ncid, "ptm_vars", &_id));
                RST_CHECK(nc_get_var_double(ncid, _id, PTM_Vars));
            }
        } else if (ptm_sw && PTM_Stat != NULL) {
            fprintf(stderr, "     WARNING: ptm_sw is on but restart has no "
                    "PTM state; particles initialised from scratch.\n");
        }
    }

    nc_close(ncid);
    rst_ncid_ = -1;
    return 1;
}
